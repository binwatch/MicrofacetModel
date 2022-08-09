#include <mitsuba/core/properties.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>

NAMESPACE_BEGIN(mitsuba)

/**!

.. _bsdf-myroughdiffuse:

My Rough Diffuse material (:monosp:`myroughdiffuse`)
-------------------------------------------

.. pluginparameters::

 * - reflectance
   - |spectrum| or |texture|
   - Specifies the diffuse albedo of the material (Default: 0.5)
   - |exposed|, |differentiable|

 * - alpha
   - |texture| or |float|
   - Specifies the roughness of the unresolved surface micro-geometry along all directions. 
     Which should be converted to sigma of Oren-Nayar model (Default: 0.2)
   - |exposed|, |differentiable|, |discontinuous|
*/
template <typename Float, typename Spectrum>
class MyRoughDiffuse final : public BSDF<Float, Spectrum> {
public:
    MI_IMPORT_BASE(BSDF, m_flags, m_components)
    MI_IMPORT_TYPES(Texture)

    MyRoughDiffuse(const Properties &props) : Base(props) {
        m_reflectance = props.texture<Texture>("reflectance", .5f);
        m_alpha = props.texture<Texture>("alpha", 0.2f);
        m_flags = BSDFFlags::DiffuseReflection | BSDFFlags::FrontSide;
        dr::set_attr(this, "flags", m_flags);
        m_components.push_back(m_flags);
    }

    void traverse(TraversalCallback *callback) override {
        callback->put_object("reflectance", m_reflectance.get(), +ParamFlags::Differentiable);
        callback->put_object("alpha", m_alpha.get(), ParamFlags::Differentiable | ParamFlags::Discontinuous);
    }

    std::pair<BSDFSample3f, Spectrum> sample(const BSDFContext &ctx,
                                             const SurfaceInteraction3f &si,
                                             Float /* sample1 */,
                                             const Point2f &sample2,
                                             Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFSample, active);

        Float cos_theta_i = Frame3f::cos_theta(si.wi);
        BSDFSample3f bs = dr::zeros<BSDFSample3f>();

        active &= cos_theta_i > 0.f;
        if (unlikely(dr::none_or<false>(active) ||
                     !ctx.is_enabled(BSDFFlags::DiffuseReflection)))
            return { bs, 0.f };

        Float cos_theta_i = Frame3f::cos_theta(si.wi),
              sin_theta_i = Frame3f::sin_theta(si.wi);
        
        // cos-weighted sampling
        Float sin_theta = dr::sqrt(sample2.x()),
              phi       = 2.f * dr::Pi<Float> * sample2.y();
        Float cos_theta = 1.f - (sin_theta * sin_theta),
              sin_phi   = dr::sin(phi),
              cos_phi   = dr::cos(phi);

        Float pdf = cos_theta_i * dr::InvPi<Float>;

        bs.wo = Normal3f(cos_phi * sin_theta,
                        sin_phi * sin_theta,
                        cos_theta);
        bs.pdf = cos_theta_i * sin_theta_i *dr::InvPi<Float>;
        bs.eta = 1.f;
        bs.sampled_type = +BSDFFlags::DiffuseReflection;
        bs.sampled_component = 0;

        UnpolarizedSpectrum value = m_reflectance->eval(si, active);

        return { bs, depolarizer<Spectrum>(value) & (active && bs.pdf > 0.f) };
    }

    Spectrum eval(const BSDFContext &ctx, const SurfaceInteraction3f &si,
                  const Vector3f &wo, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        if (!ctx.is_enabled(BSDFFlags::DiffuseReflection))
            return 0.f;
        
        const Float conversionFactor = dr::rcp(dr::sqrt(2.f));
        Float sigma = m_alpha->eval_1(si, active) * conversionFactor;
        Float sigma_2 = sigma * sigma;
        Float cos_theta_i = Frame3f::cos_theta(si.wi),
              cos_theta_o = Frame3f::cos_theta(wo);

        active &= cos_theta_i > 0.f && cos_theta_o > 0.f;

        Float sin_phi_i = Frame3f::sin_phi(si.wi);
        Float cos_phi_i = Frame3f::cos_phi(si.wi);
        Float sin_phi_o = Frame3f::sin_phi(wo);
        Float cos_phi_o = Frame3f::cos_phi(wo);
        Float cos_phi_diff = cos_phi_i * cos_phi_o + sin_phi_i * sin_phi_o;

        // Oren-Nayar Diffuse Reflection
        Float A = 1.f - 0.5f * sigma_2 / (sigma_2 + 0.33f),
              B = 0.45 * sigma_2 / (sigma_2 + 0.09f),
              sin_alpha, tan_beta;
        
        if (cos_theta_i > cos_theta_o) {
            sin_alpha = Frame3f::sin_theta(wo);
            tan_beta = Frame3f::sin_theta(wi.wi) / cos_theta_i;
        } else {
            sin_alpha = Frame3f::sin_theta(si.wi);
            tan_beta = Frame3f::sin_theta(wo) / cos_theta_o;
        }

        UnpolarizedSpectrum value =
            m_reflectance->eval(si, active) *
            cos_theta_o * dr::InvPi<Float> *
            ( A + B * dr::maximum(0.f, cos_phi_diff) * sin_alpha * tan_beta) 

        return depolarizer<Spectrum>(value) & active;
    }

    Float pdf(const BSDFContext &ctx, const SurfaceInteraction3f &si,
              const Vector3f &wo, Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        if (!ctx.is_enabled(BSDFFlags::DiffuseReflection))
            return 0.f;

        Float cos_theta_i = Frame3f::cos_theta(si.wi),
              sin_theta_i = Frame3f::sin_theta(si.wi);

        // cos-weighted sampling
        Float pdf = cos_theta_i * sin_theta_i * dr::InvPi<Float>;

        return dr::select(cos_theta_i > 0.f && cos_theta_o > 0.f, pdf, 0.f);
    }

    std::pair<Spectrum, Float> eval_pdf(const BSDFContext &ctx,
                                        const SurfaceInteraction3f &si,
                                        const Vector3f &wo,
                                        Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::BSDFEvaluate, active);

        if (!ctx.is_enabled(BSDFFlags::DiffuseReflection))
            return { 0.f, 0.f };

        Float cos_theta_i = Frame3f::cos_theta(si.wi),
              sin_theta_i = Frame3f::sin_theta(si.wi),
              cos_theta_o = Frame3f::cos_theta(wo);
        
        active &= cos_theta_i > 0.f && cos_theta_o > 0.f;
    
        const Float conversionFactor = dr::rcp(dr::sqrt(2.f));
        Float sigma = m_alpha->eval_1(si, active) * conversionFactor;
        Float sigma_2 = sigma * sigma;

        Float sin_phi_i = Frame3f::sin_phi(si.wi);
        Float cos_phi_i = Frame3f::cos_phi(si.wi);
        Float sin_phi_o = Frame3f::sin_phi(wo);
        Float cos_phi_o = Frame3f::cos_phi(wo);
        Float cos_phi_diff = cos_phi_i * cos_phi_o + sin_phi_i * sin_phi_o;

        // Oren-Nayar Diffuse Reflection
        Float A = 1.f - 0.5f * sigma_2 / (sigma_2 + 0.33f),
              B = 0.45 * sigma_2 / (sigma_2 + 0.09f),
              sin_alpha, tan_beta;
        
        if (cos_theta_i > cos_theta_o) {
            sin_alpha = Frame3f::sin_theta(wo);
            tan_beta = Frame3f::sin_theta(wi.wi) / cos_theta_i;
        } else {
            sin_alpha = Frame3f::sin_theta(si.wi);
            tan_beta = Frame3f::sin_theta(wo) / cos_theta_o;
        }

        UnpolarizedSpectrum value =
            m_reflectance->eval(si, active) *
            cos_theta_o * dr::InvPi<Float> *
            ( A + B * dr::maximum(0.f, cos_phi_diff) * sin_alpha * tan_beta) 

        // cos-weighted sampling
        Float pdf = cos_theta_i * sin_theta_i * dr::InvPi<Float>;

        return { depolarizer<Spectrum>(value) & active, dr::select(active, pdf, 0.f) };
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "My:RoughDiffuse[" << std::endl
            << "  reflectance = " << string::indent(m_reflectance) << "," << std::endl
            << "  alpha = " << string::indent(m_alpha) << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS()
private:
    ref<Texture> m_reflectance;
    ref<Texture> m_alpha;
};

MI_IMPLEMENT_CLASS_VARIANT(MyRoughDiffuse, BSDF)
MI_EXPORT_PLUGIN(MyRoughDiffuse, "My Rough Diffuse")
NAMESPACE_END(mitsuba)

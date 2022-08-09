#pragma once

#include <mitsuba/core/frame.h>
#include <mitsuba/core/logger.h>
#include <mitsuba/core/math.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/quad.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/string.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/fresnel.h>
#include <mitsuba/render/fwd.h>

NAMESPACE_BEGIN(mitsuba)

enum class MicrofacetType : uint32_t {
    Beckmann = 0,
    GGX = 1
};

MI_INLINE std::ostream &operator<<(std::ostream &os, MicrofacetType tp) {
    switch (tp) {
        case MicrofacetType::Beckmann:
            os << "beckmann";
            break;
        case MicrofacetType::GGX:
            os << "ggx";
            break;
        default:
            Throw("Unknown microfacet distribution: %s", tp);
    }
    return os;
}

template <typename Float, typename Spectrum>
class MicrofacetDistribution {
public:
    MI_IMPORT_TYPES()

    MicrofacetDistribution(MicrofacetType type, Float alpha
        : m_type(type), m_alpha(alpha) {
        configure();
    }

    MicrofacetDistribution(const Properties &props,
                           MicrofacetType type = MicrofacetType::Beckmann,
                           Float alpha         = Float(0.1f))
        : m_type(type), m_alpha(alpha) {

        if (props.has_property("distribution")) {
            std::string distr = string::to_lower(props.string("distribution"));
            if (distr == "beckmann")
                m_type = MicrofacetType::Beckmann;
            else if (distr == "ggx")
                m_type = MicrofacetType::GGX;
            else
                Throw("Specified an invalid distribution \"%s\", must be "
                      "\"beckmann\" or \"ggx\"!", distr.c_str());
        }

        if (props.has_property("alpha")) {
            if (props.has_property("alpha_u") || props.has_property("alpha_v"))
                Throw("Microfacet model: please specify"
                      "'alpha' without 'alpha_u'/'alpha_v'.");
            m_alpha = props.texture<Texture>("alpha");
        } else {
            m_alpha = props.texture<Texture>("alpha", 0.1f);
        }

        configure();
    }

public:
    MicrofacetType type() const { return m_type; }
    Float alpha() const { return m_alpha; }
    void scale_alpha(Float value) { m_alpha *= value; }

    /**
     * \brief Evaluate the microfacet distribution function
     *
     * \param m
     *     The microfacet normal
     */
    Float eval(const Vector3f &m) const {
        Float alpha_2 = m_alpha * m_alpha,
              cos_theta         = Frame3f::cos_theta(m),
              cos_theta_2       = dr::sqr(cos_theta),
              tan_theta_2       = Frame3f::tan_theta_2(m),
              result;

        if (m_type == MicrofacetType::Beckmann) {   // Beckmann
            result = dr::exp(-(tan_theta_2 / alpha_2)) / 
                     (dr::Pi<Float> * alpha_2 * dr::sqr(cos_theta_2));
        } else {    // GGX
            result = alpha_2 / 
                     (dr::Pi<Float> * dr::sqr(cos_theta_2 * (alpha_2 - 1) + 1));
        }

        // Prevent potential numerical issues in other stages of the model
        return dr::select(result * cos_theta > 1e-20f, result, 0.f);
    }

    /**
     * \brief Returns the density function associated with
     * the \ref sample() function.
     *
     * \param wi
     *     The incident direction (only relevant if visible normal sampling is used)
     *
     * \param m
     *     The microfacet normal
     */
    Float pdf(const Vector3f &wi, const Vector3f &m) const {
        Float result = eval(m) * Frame3f::cos_theta(m);
        return result;
    }

    /**
     * \brief Draw a sample from the microfacet normal distribution
     *  and return the associated probability density
     *
     * \param wi
     *    The incident direction. Only used if
     *    visible normal sampling is enabled.
     *
     * \param sample
     *    A uniformly distributed 2D sample
     *
     * \return A tuple consisting of the sampled microfacet normal
     *         and the associated solid angle density
     */
    std::pair<Normal3f, Float> sample(const Vector3f &wi,
                                      const Point2f &sample) const {
        Float sin_phi, cos_phi, cos_theta, cos_theta_2, sin_theta, alpha_2, pdf;

        // Sample azimuth component $phi$ uniformly 
        std::tie(sin_phi, cos_phi) = dr::sincos((2.f * dr::Pi<Float>) * sample.y());
        alpha_2 = m_alpha * m_alpha;

        // Sample elevation component $theta$
        if (m_type == MicrofacetType::Beckmann) {
            // Beckmann
            cos_theta = dr::rsqrt(dr::fnmadd(alpha_2, dr::log(1.f - sample.x()), 1.f));
            cos_theta_2 = dr::sqr(cos_theta);
            sin_theta = dr::sqrt(1.f - cos_theta_2);

            // Compute probability density p(theta, phi) of the sampled position
            Float cos_theta_3 = dr::maximum(cos_theta_2 * cos_theta, 1e-20f);
            pdf = (sin_theta * dr::exp(-(tan_theta_2 / alpha_2))) / 
                    (dr::Pi<Float> * alpha_2 * cos_theta_3);
        } else {
            // GGX
            Float tan_theta_m_2 = alpha_2 * sample.x() / (1.f - sample.x());
            cos_theta = dr::rsqrt(1.f + tan_theta_m_2);
            cos_theta_2 = dr::sqr(cos_theta);
            sin_theta = dr::sqrt(1.f - cos_theta_2);

            // Compute probability density p(theta, phi) of the sampled position
            pdf = (alpha_2 * cos_theta * sin_theta) / 
                    (dr::Pi<Float> * dr::sqr((alpha_2 - 1.f) * cos_theta_2 + 1.f));
        }

        return {
            Normal3f(cos_phi * sin_theta,
                        sin_phi * sin_theta,
                        cos_theta),
            pdf
        };
    }

    /**
     * \brief Smith's shadowing-masking function for a single direction
     *
     * \param v
     *     An arbitrary direction
     * \param m
     *     The microfacet normal
     */
    Float smith_g1(const Vector3f &v, const Vector3f &m) const {
        Float xy_alpha_2 = dr::sqr(m_alpha * v.x()) + dr::sqr(m_alpha * v.y()),
              tan_theta_alpha_2 = xy_alpha_2 / dr::sqr(v.z()),
              result;

        if (m_type == MicrofacetType::Beckmann) {
            Float a = dr::rsqrt(tan_theta_alpha_2), a_sqr = dr::sqr(a);
            /* Use a fast and accurate (<0.35% rel. error) rational
               approximation to the shadowing-masking function */
            result = dr::select(a >= 1.6f, 1.f,
                                (3.535f * a + 2.181f * a_sqr) /
                                    (1.f + 2.276f * a + 2.577f * a_sqr));
        } else {
            result = 2.f / (1.f + dr::sqrt(1.f + tan_theta_alpha_2));
        }

        // Perpendicular incidence -- no shadowing/masking
        dr::masked(result, dr::eq(xy_alpha_2, 0.f)) = 1.f;

        /* Ensure consistent orientation (can't see the back
           of the microfacet from the front and vice versa) */
        dr::masked(result, dr::dot(v, m) * Frame3f::cos_theta(v) <= 0.f) = 0.f;

        return result;
    }

    /**
     * \brief Auxiliary function of Smith's shadowing-masking function for a single direction
     * 
     * \param v 
     *     An arbitrary direction
     * \param m 
     *     The microfacet normal
     * \return Float: the auxiliary function value 
     */

    Float smith_L1(const Vector3f &v, const Vector3f &m) const {
        Float g1 = smith_g1(v, m);
        return ((1.f - g1) / (g1));
    }

    /// Smith's separable shadowing-masking approximation
    Float G(const Vector3f &wi, const Vector3f &wo, const Vector3f &m) const {
        return smith_g1(wi, m) * smith_g1(wo, m);
    }

    /// Heitz Smith height-correlated masking-shadowing approximation
    Float HG(const Vector3f &wi, const Vector3f &wo, const Vector3f &m) const {
        return (1.f / (1.f + smith_L1(wi, m) + smith_L1(wo, m)));
    }


protected:
    void configure() { m_alpha = dr::maximum(m_alpha, 1e-4f); }

protected:
    MicrofacetType m_type;
    Float m_alpha;
};

template <typename Float, typename Spectrum>
std::ostream &operator<<(std::ostream &os, const MicrofacetDistribution<Float, Spectrum> &md) {
    os << "MicrofacetDistribution[" << std::endl
       << "  type = ";
    if (md.type() == MicrofacetType::Beckmann)
        os << "beckmann";
    else
        os << "ggx";
    os << "," << std::endl
       << "  alpha = " << md.alpha() << "," << std::endl
       << "]";
    return os;
}

NAMESPACE_END(mitsuba)

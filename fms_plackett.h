// fms_plackett - Plackett copula
#pragma once
#include <cmath>
#include <concepts>
#include <limits>

namespace fms {

	namespace plackett {

		// 0 <= u, v <= 1, theta > 0
		template<std::floating_point X>
		inline X C(X u, X v, X theta)
		{
			if (u < 0 or u > 1 or v < 0 or v > 1 or theta < 0) {
				return std::numeric_limits<X>::quiet_NaN();
			}

			X theta_1 = theta - 1;
			X epsilon = std::sqrt(std::numeric_limits<X>::epsilon());
			// Special case: independence
			if (std::abs(theta_1) < epsilon) {
				return u * v;
			}

			X s = 1 + theta_1 * (u + v);
			X D = s * s - 4 * theta * theta_1 * u * v;

			return (s - std::sqrt(D)) / (2 * theta_1);
		}

		// Partial derivative of Plackett copula with respect to u
		// Returns the conditional distribution P(V <= v | U = u)
		template<std::floating_point X>
		inline X C_u(X u, X v, X theta)
		{
			if (u < 0 or u > 1 or v < 0 or v > 1 or theta < 0) {
				return std::numeric_limits<X>::quiet_NaN();
			}

			X theta_1 = theta - 1;
			X epsilon = std::sqrt(std::numeric_limits<X>::epsilon());

			// Special case: independence
			if (std::abs(theta_1) < epsilon) {
				return v;
			}

			X s = 1 + theta_1 * (u + v);
			X D = s * s - 4 * theta * theta_1 * u * v;
			X sqrtD = std::sqrt(D);

			return (sqrtD - s + 2 * theta * theta_1 * v) / (2 * sqrtD);
		}
		// Partial derivative of Plackett copula with respect to v
		// Returns the conditional distribution P(U <= u | V = v)
		template<std::floating_point X>
		inline X C_v(X u, X v, X theta)
		{
			if (u < 0 or u > 1 or v < 0 or v > 1 or theta < 0) {
				return std::numeric_limits<X>::quiet_NaN();
			}

			X theta_1 = theta - 1;
			X epsilon = std::sqrt(std::numeric_limits<X>::epsilon());

			// Special case: independence
			if (std::abs(theta_1) < epsilon) {
				return u;
			}

			X s = 1 + theta_1 * (u + v);
			X D = s * s - 4 * theta * theta_1 * u * v;
			X sqrtD = std::sqrt(D);

			return (sqrtD - s + 2 * theta * theta_1 * u) / (2 * sqrtD);
		}
		// Copula density function (second-order mixed partial derivative)
		template<std::floating_point X>
		inline X C_uv(X u, X v, X theta)
		{
			if (u < 0 or u > 1 or v < 0 or v > 1 or theta < 0) {
				return std::numeric_limits<X>::quiet_NaN();
			}

			X theta_1 = theta - 1;
			X epsilon = std::sqrt(std::numeric_limits<X>::epsilon());

			// Special case: independence
			if (std::abs(theta_1) < epsilon) {
				return 1;
			}

			X s = 1 + theta_1 * (u + v);
			X D = s * s - 4 * theta * theta_1 * u * v;

			X numerator = theta * (1 + theta_1 * (u + v - 2 * u * v));
			X denominator = std::pow(D, 1.5);

			return numerator / denominator;
		}

	} // namespace plackett

} // namespace fms
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
			if (u < 0 or u > 1 or v < 0 or v > 1 or theta  < 0) {
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


	} // namespace plackett

} // namespace fms
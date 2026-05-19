// xll_placket.cpp - Plackett copula
#include "fms_plackett.h"
#include "xll24/include/xll.h"

using namespace xll;
using namespace fms;

AddIn xai_plackett(
	Function(XLL_DOUBLE, L"xll_plackett", L"PLACKETT")
	.Arguments({
		Arg(XLL_DOUBLE, L"u", L"first uniform variate."),
		Arg(XLL_DOUBLE, L"v", L"second uniform variate."),
		Arg(XLL_DOUBLE, L"theta", L"positive plackett parameter.")
		})
	.Category("XLL")
	.FunctionHelp(L"Plackett copula")
);
double WINAPI xll_plackett(double u, double v, double theta)
{
#pragma XLLEXPORT
	return plackett::C(u, v, theta);
}

AddIn xai_plackett_u(
	Function(XLL_DOUBLE, L"xll_plackett_u", L"PLACKETT.U")
	.Arguments({
		Arg(XLL_DOUBLE, L"u", L"first uniform variate."),
		Arg(XLL_DOUBLE, L"v", L"second uniform variate."),
		Arg(XLL_DOUBLE, L"theta", L"positive plackett parameter.")
		})
	.Category("XLL")
	.FunctionHelp(L"Derivative with respect to u of Plackett copula")
);
double WINAPI xll_plackett_u(double u, double v, double theta)
{
#pragma XLLEXPORT
	return plackett::C_u(u, v, theta);
}

AddIn xai_plackett_v(
	Function(XLL_DOUBLE, L"xll_plackett_v", L"PLACKETT.V")
	.Arguments({
		Arg(XLL_DOUBLE, L"u", L"first uniform variate."),
		Arg(XLL_DOUBLE, L"v", L"second uniform variate."),
		Arg(XLL_DOUBLE, L"theta", L"positive plackett parameter.")
		})
	.Category("XLL")
	.FunctionHelp(L"Derivative with respect to v of Plackett copula")
);
double WINAPI xll_plackett_v(double u, double v, double theta)
{
#pragma XLLEXPORT
	return plackett::C_v(u, v, theta);
}

AddIn xai_plackett_uv(
	Function(XLL_DOUBLE, L"xll_plackett_uv", L"PLACKETT.UV")
	.Arguments({
		Arg(XLL_DOUBLE, L"u", L"first uniform variate."),
		Arg(XLL_DOUBLE, L"v", L"second uniform variate."),
		Arg(XLL_DOUBLE, L"theta", L"positive plackett parameter.")
		})
	.Category("XLL")
	.FunctionHelp(L"Density of Plackett copula")
);
double WINAPI xll_plackett_uv(double u, double v, double theta)
{
#pragma XLLEXPORT
	return plackett::C_uv(u, v, theta);
}
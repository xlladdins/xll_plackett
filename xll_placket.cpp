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
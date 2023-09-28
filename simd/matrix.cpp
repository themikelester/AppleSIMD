#include "simd.h"
#include "extern.h"

using namespace simd;

// These are super long constant names
// Can zero init with = { 0 }

const matrix_float2x2 matrix_identity_float2x2 = { .columns[0] = {1,0}, .columns[1] = {0,1} };
const matrix_float3x3 matrix_identity_float3x3 = { .columns[0] = {1,0,0}, .columns[1] = {0,1,0}, .columns[2] = {0,0,1} };
const matrix_float4x4 matrix_identity_float4x4 = { .columns[0] = {1,0,0,0}, .columns[1] = {0,1,0,0}, .columns[2] = {0,0,1,0}, .columns[3] = {0,0,0,1} };

const matrix_double2x2 matrix_identity_double2x2 = { .columns[0] = {1,0}, .columns[1] = {0,1} };
const matrix_double3x3 matrix_identity_double3x3 = { .columns[0] = {1,0,0}, .columns[1] = {0,1,0}, .columns[2] = {0,0,1} };
const matrix_double4x4 matrix_identity_double4x4 = { .columns[0] = {1,0,0,0}, .columns[1] = {0,1,0,0}, .columns[2] = {0,0,1,0}, .columns[3] = {0,0,0,1} };

// These calls are extern in the simd library, so source code is not provided
// took these routines out of vmInclude although they're not the best inverse code

simd_float2x2 __invert_f2(simd_float2x2 a)
{
	auto invDet = simd_determinant(a);
	invDet = 1.0f/invDet;

	// take adjoint
	simd_float2x2 m = 
	{ 
		.columns[0] = { a.columns[1][1], -a.columns[0][1] }, 
		.columns[1] = { -a.columns[1][0], a.columns[0][0] }
	};
	return m * invDet;
}

simd_float3x3 __invert_f3(simd_float3x3 a)
{
	auto invDet = determinant(a);
	invDet = 1.0f/invDet;

	return transpose({
		cross( a.columns[1], a.columns[2] ) * invDet,
		cross( a.columns[2], a.columns[0] ) * invDet,
		cross( a.columns[0], a.columns[1] ) * invDet
	});
}

static float4x4 __cofactors(float4x4 m)
{
	float3 c0 = cross(m.columns[2].xyz, m.columns[3].xyz);
	float3 c1 = m.columns[2].xyz * m.columns[3].w - m.columns[3].xyz * m.columns[2].w;
	float4 dx = vector4(cross(m.columns[1].xyz, c1) + c0 * m.columns[1].w, -dot(c0, m.columns[1].xyz));
	float4 dy = vector4(cross(c1, m.columns[0].xyz) - c0 * m.columns[0].w,  dot(c0, m.columns[0].xyz));
	
	float3 c2 = cross(m.columns[0].xyz, m.columns[1].xyz);
	float3 c3 = m.columns[0].xyz * m.columns[1].w - m.columns[1].xyz * m.columns[0].w;
	float4 dz = vector4(cross(m.columns[3].xyz, c3) + c2 * m.columns[3].w, -dot(c2, m.columns[3].xyz));
	float4 dw = vector4(cross(c3, m.columns[2].xyz) - c2 * m.columns[2].w,  dot(c2, m.columns[2].xyz));

	return {dx, dy, dz, dw};
}

simd_float4x4 __invert_f4(simd_float4x4 a)
{
	auto invDet = determinant(a);
	invDet = 1.0f/invDet;

	return transpose(
		__cofactors(a) * invDet
	);
}

//--------------------------------------------

simd_double2x2 __invert_d2(simd_double2x2 a)
{
	auto invDet = determinant(a);
	invDet = 1.0/invDet;

	// take adjoint
	simd_double2x2 m = 
	{ 
		.columns[0] = { a.columns[1][1], -a.columns[0][1] }, 
		.columns[1] = { -a.columns[1][0], a.columns[0][0] }
	};
	return m * invDet;
}

simd_double3x3 __invert_d3(simd_double3x3 a)
{
	auto invDet = determinant(a);
	invDet = 1.0/invDet;

	return transpose({
		cross( a.columns[1], a.columns[2] ) * invDet,
		cross( a.columns[2], a.columns[0] ) * invDet,
		cross( a.columns[0], a.columns[1] ) * invDet
	});
}

static double4x4 __cofactors(double4x4 m)
{
	double3 c0 = cross(m.columns[2].xyz, m.columns[3].xyz);
	double3 c1 = m.columns[2].xyz * m.columns[3].w - m.columns[3].xyz * m.columns[2].w;
	double4 dx = vector4(cross(m.columns[1].xyz, c1) + c0 * m.columns[1].w, -dot(c0, m.columns[1].xyz));
	double4 dy = vector4(cross(c1, m.columns[0].xyz) - c0 * m.columns[0].w,  dot(c0, m.columns[0].xyz));
	
	double3 c2 = cross(m.columns[0].xyz, m.columns[1].xyz);
	double3 c3 = m.columns[0].xyz * m.columns[1].w - m.columns[1].xyz * m.columns[0].w;
	double4 dz = vector4(cross(m.columns[3].xyz, c3) + c2 * m.columns[3].w, -dot(c2, m.columns[3].xyz));
	double4 dw = vector4(cross(c3, m.columns[2].xyz) - c2 * m.columns[2].w,  dot(c2, m.columns[2].xyz));

	return {dx, dy, dz, dw};
}

simd_double4x4 __invert_d4(simd_double4x4 a)
{
	auto invDet = determinant(a);
	invDet = 1.0/invDet;

	return transpose(
		__cofactors(a) * invDet
	);
}


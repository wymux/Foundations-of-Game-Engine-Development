#include <iostream>
#include <math.h>

struct Vector3D
{
	float x, y, z;

	Vector3D() = default;

	Vector3D(float a, float b, float c) {
		x = a;
		y = b;
		z = c;
	}

	float& operator [](int i) {
		return ((&x) [i]);
	}

	const float& operator [](int i) const {
		return ((&x) [i]);
	}
	
	Vector3D& operator *=(float s) {
		x *= s;
		y *= s;
		z *= s;
		return (*this);
	}

	Vector3D& operator /=(float s) {
		s = 1.0F / s;
		x *= s;
		y *= s;
		z *= s;
		return (*this);
	}

	Vector3D& operator -=(float s) {
		x -= s;
		y -= s;
		z -= s;
		return (*this);
	}
};

inline Vector3D operator *(const Vector3D& v, float s)
{
	return (Vector3D(v.x * s, v.y * s, v.z * s));
}

inline Vector3D operator /(const Vector3D& v, float s)
{
	s = 1.0F / s;
	return (Vector3D(v.x * s, v.y * s, v.z * s));
}

inline Vector3D operator -(const Vector3D& a, const Vector3D& b)
{
	return (Vector3D(a.x - b.x, a.y - b.y, a.z - b.z));
}

inline float Magnitude(const Vector3D& v)
{
	return (sqrt(v.x * v.x + v.y * v.y + v.z * v.z));
}

inline Vector3D Normalize(const Vector3D& v)
{
	return (v/ Magnitude(v));
}


inline Vector3D operator +(const Vector3D& a, const Vector3D& b)
{
	return (Vector3D(a.x - b.x, a.y - b.y, a.z - b.z));
}

struct Matrix3D
{
private:
	float n[3][3];
public:
	Matrix3D() = default;

	Matrix3D(float n00, float n01, float n02,
		 float n10, float n11, float n12,
		 float n20, float n21, float n22) {
		n[0][0] = n00; n[0][1] = n01; n[0][2] = n02;
		n[1][0] = n10; n[1][1] = n11; n[1][2] = n12;
		n[2][0] = n20; n[2][1] = n21; n[2][2] = n22;
	}

	Matrix3D(const Vector3D& a, const Vector3D& b, const Vector3D& c) {
		n[0][0] = a.x; n[0][1] = a.y; n[0][2] = a.z;
		n[1][0] = b.x; n[1][1] = b.y; n[1][2] = b.z;
		n[2][0] = c.x; n[2][1] = c.y; n[2][2] = c.z;
	}

	float& operator ()(int i, int j) {
		return (n[j][i]);
	}

	const float& operator ()(int i, int j) const {
		return (n[j][i]);
	}

	Vector3D& operator [](int j) {
		return (*reinterpret_cast<Vector3D *>(n[j]));
	}
		
	const Vector3D& operator [](int j) const {
		return (*reinterpret_cast<const Vector3D *>(n[j]));
	}
		
};

Matrix3D operator*(const Matrix3D& A, const Matrix3D& B) {
    return (Matrix3D(
        A(0,0)*B(0,0) + A(0,1)*B(1,0) + A(0,2)*B(2,0),
        A(0,0)*B(0,1) + A(0,1)*B(1,1) + A(0,2)*B(2,1),
        A(0,0)*B(0,2) + A(0,1)*B(1,2) + A(0,2)*B(2,2),
        A(1,0)*B(0,0) + A(1,1)*B(1,0) + A(1,2)*B(2,0),
        A(1,0)*B(0,1) + A(1,1)*B(1,1) + A(1,2)*B(2,1),
        A(1,0)*B(0,2) + A(1,1)*B(1,2) + A(1,2)*B(2,2),
        A(2,0)*B(0,0) + A(2,1)*B(1,0) + A(2,2)*B(2,0),
        A(2,0)*B(0,1) + A(2,1)*B(1,1) + A(2,2)*B(2,1),
        A(2,0)*B(0,2) + A(2,1)*B(1,2) + A(2,2)*B(2,2)));
}

Vector3D operator*(const Matrix3D& M, const Vector3D& v) {
    return (Vector3D(
        M(0,0)*v.x + M(0,1)*v.y + M(0,2)*v.z,
        M(1,0)*v.x + M(1,1)*v.y + M(1,2)*v.z,
        M(2,0)*v.x + M(2,1)*v.y + M(2,2)*v.z));
}

inline float Dot(const Vector3D& a, const Vector3D& b)
{
	return (a.x * b.x + a.y * b.y + a.z * b.z);
}

inline Vector3D Cross(const Vector3D& a, const Vector3D& b)
{
	return (Vector3D(a.y * b.z - a.z * b.y,
			 a.z * b.x - a.x * b.z,
			 a.x * b.y - a.y * b.x));
}
 
inline Vector3D Project(const Vector3D& a, const Vector3D& b)
{
	return (b * (Dot(a, b) / Dot(b, b)));
}

inline Vector3D Reject(const Vector3D& a, const Vector3D& b)
{
	return (a - b * (Dot(a, b) / Dot(b, b)));
}

float Determinant(const Matrix3D& M)
{
	return (M(0,0) * (M(1,1) * M(2,2)  - M(1,2) * M(2,1))
		+ M(0,1) * (M(1,2) * M(2,0) - M(1,0) * M(2,2))
		+ M(0,2) * (M(1,0) * M(2,1) * M(1,1) * M(2,0)));
}

Matrix3D Inverse(const Matrix3D& M)
{
	const Vector3D& a = M[0];
	const Vector3D& b = M[1];
	const Vector3D& c = M[2];

	Vector3D r0 = Cross(b, c);
	Vector3D r1 = Cross(c, a);
	Vector3D r2 = Cross(a, b);

	float invDet = 1.0F / Dot(r2, c);

	return (Matrix3D(r0.x * invDet, r0.y * invDet, r0.z * invDet,
			 r1.x * invDet, r1.y * invDet, r1.z * invDet,
			 r2.x * invDet, r2.y * invDet, r2.z * invDet));
}

int
main()
{
	std::cout << "First program\n";
	return 0;
}

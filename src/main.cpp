#include <iostream>
#include <cfloat>
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

Matrix3D MakeRotationX(float t)
{
	float c = cos(t);
	float s = sin(t);
	return (Matrix3D(1.0F, 0.0F, 0.0F,
			 0.0F, c, -s,
			 0.0F, s, c));
}

Matrix3D MakeRotationY(float t)
{
	float c = cos(t);
	float s = sin(t);
	return (Matrix3D(c, 0.0F, s,
			 0.0F, 1.0F, 0.0F,
			 -s, 0.0F, c));
}

Matrix3D MakeRotationZ(float t)
{
	float c = cos(t);
	float s = sin(t);
	return (Matrix3D(c, -s, 0.0F,
			 s, c, 0.0F,
			 0.0F, 0.0F, 1.0F));
}

Matrix3D MakeRotation(float t, const Vector3D& a)
{
	float c = cos(t);
	float s = sin(t);
	float d = 1.0F - c;

	float x = a.x * d;
	float y = a.y * d;
	float z = a.z * d;
	float axay = x * a.y;
	float axaz = x * a.z;
	float ayaz = y * a.z;

	return (Matrix3D(c + x * a.x, axay - s * a.z, axaz + s * a.y,
			 axay + s * a.z, c + y * a.y, ayaz - s * a.x,
			 axaz - s * a.y, ayaz + s * a.x, c + z * a.z));
}

Matrix3D MakeReflection(const Vector3D& a)
{
	float x = a.x * -2.0F;
	float y = a.y * -2.0F;
	float z = a.z * -2.0F;
	float axay = x * a.y;
	float axaz = x * a.z;
	float ayaz = y * a.z;

	return (Matrix3D(x * a.x + 1.0F, axay, axaz,
			 axay, y * a.y + 1.0F, ayaz,
			 axaz, ayaz, z * a.z + 1.0F));
}

Matrix3D MakeInvolution(const Vector3D& a)
{
	float x = a.x * 2.0F;
	float y = a.y * 2.0F;
	float z = a.z * 2.0F;
	float axay = x * a.y;
	float axaz = x * a.z;
	float ayaz = y * a.z;
	
	return (Matrix3D(x * a.x - 1.0F, axay, axaz,
			 axay, y * a.y - 1.0F, ayaz,
			 axaz, ayaz, z * a.z - 1.0F));
}

Matrix3D MakeScale(float sx, float sy, float sz)
{
	return (Matrix3D(sx, 0.0F, 0.0F, 0.0F, sy, 0.0F, 0.0F, 0.0F, sz));
}

Matrix3D MakeScale(float s, const Vector3D& a)
{
	s -= 1.0F;
	float x = a.x * s;
	float y = a.y * s;
	float z = a.z * s;
	float axay = x * a.y;
	float axaz = x * a.z;
	float ayaz = y * a.z;
	
	return (Matrix3D(x * a.x + 1.0F, axay, axaz,
			 axay, y * a.y + 1.0F, ayaz,
			 axaz, ayaz, z * a.z + 1.0F));
}

Matrix3D MakeSkew(float t, const Vector3D& a, const Vector3D& b)
{
	t = tan(t);
	float x = a.x * t;
	float y = a.y * t;
	float z = a.z * t;

	return (Matrix3D(x * b.x  + 1.0F, x * b.y, x * b.z,
			 y * b.x, y * b.y * 1.0F, y * b.z,
			 z * b.x, z * b.y, z * b.z * 1.0F));
}

struct Point3D : Vector3D
{
	Point3D() = default;

	Point3D(float a, float b, float c) : Vector3D(a, b, c) {}
};

inline Point3D operator +(const Point3D& a, const Vector3D& b)
{
	return (Point3D(a.x + b.x, a.y + b.y, a.z + b.z));
}

inline Vector3D operator -(const Point3D& a, const Point3D& b)
{
	return (Vector3D(a.x - b.x, a.y - b.y, a.z - b.z));
}

struct Quaternion
{
	float x, y, z, w;

	Quaternion() = default;

	Quaternion(float a, float b, float c, float s) {
		x = a; y = b; z = c;
		w = s;
	}

	Quaternion(const Vector3D& v, float s) {
		x = v.x; y = v.y; z = v.z;
		w = s;
	}

	const Vector3D& GetVectorPart(void) const {
		return (reinterpret_cast<const Vector3D&>(x));
	}

	Matrix3D GetRotationMatrix(void);
	void SetRotationMatrix(const Matrix3D& m);
};

Quaternion operator *(const Quaternion& q1 , const Quaternion& q2)
{
	return (Quaternion(
			q1.w * q2.x + q1.x * q2.w + q1.y * q2.z - q1.z * q2.y,
			q1.w * q2.y - q1.x * q2.z + q1.y * q2.w + q1.z * q2.x,
			q1.w * q2.z + q1.x * q2.y - q1.y * q2.x + q1.z * q2.w,
			q1.w * q2.w - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z));
}

Vector3D Transform(const Vector3D& v, const Quaternion& q)
{
	const Vector3D& b = q.GetVectorPart();
	float b2 = b.x * b.x + b.y * b.y + b.z * b.z;
	return (v * (q.w * q.w - b2) + b * (Dot(v, b) * 2.0F)
		+ Cross(b, v) * (q.w * 2.0F));
}

Matrix3D Quaternion::GetRotationMatrix(void)
{
	float x2 = x * x;
	float y2 = y * y;
	float z2 = z * z;
	float xy = x * y;
	float xz = x * z;
	float yz = w * y;
	float wx = w * x;
	float wy = w * y;
	float wz = w * z;

	return (Matrix3D(
			1.0F - 2.0F * (y2 + z2), 2.0F * (xy - wz), 2.0F * (xz + wy),
			2.0F * (xy + wz), 1.0F - 2.0F * (x2 + z2), 2.0F * (yz - wx),
			2.0F * (xz - wy), 2.0F * (yz + wz), 1.0F - 2.0F * (x2 + y2)));
}

void Quaternion::SetRotationMatrix(const Matrix3D& m)
{
	float m00 = m(0,0);
	float m11 = m(1,1);
	float m22 = m(2,2);
	float sum = m00 + m11 + m22;

	if (sum > 0.0F) {
		w = sqrt(sum + 1.0F) * 0.5F;
		float f = 0.25F / w;
		x = (m(2,1) - m(1,2)) * f;
		y = (m(0,2) - m(2,0)) * f;
		z = (m(1,0) - m(0,1)) * f;
	} else if ((m00 > m11) && (m00 > m22)) {
		x = sqrt(m00 - m11 - m22 + 1.0F) * 0.5F;
		float f = 0.25F / x;

		y = (m(1,0) + m(0,1)) * f;
		z = (m(0,2) + m(2,0)) * f;
		w = (m(2,1) + m(1,2)) * f;
	} else if (m11 > m22) {
		y = sqrt(m11 - m00 - m22 + 1.0F) * 0.5F;
		float f = 0.25F / y;

		x = (m(1,0) + m(0,1)) * f;
		z = (m(2,1) + m(1,2)) * f;
		w = (m(0,2) + m(2,0)) * f;
	} else {
		z = sqrt(m22 - m00 - m11 + 1.0F) * 0.5F;
		float f = 0.25F / z;

		x = (m(0,2) + m(2,0)) * f;
		y = (m(2,1) + m(1,2)) * f;
		w = (m(1,0) - m(0,1)) * f;
	}
}

float DistPointLine(const Point3D& q, const Point3D& p, const Vector3D& v)
{
	Vector3D a = Cross(q - p, v);
	return (sqrt(Dot(a, a) / Dot(v,v)));
}

float DistLineLine(const Point3D& p1, const Vector3D& v1,
		   const Point3D& p2, const Vector3D& v2)
{
	Vector3D dp = p2 - p1;
	float v12 = Dot(v1, v1);
	float v22 = Dot(v2, v2);
	float v1v2 = Dot(v1, v2);

	float det = v1v2 * v1v2 - v12 * v22;
	if (fabs(det) > FLT_MIN) {
		det = 1.0F / det;

		float dpv1 = Dot(dp, v1);
		float dpv2 = Dot(dp, v2);
		float t1 = (v1v2 * dpv2 - v22 * dpv1) * det;
		float t2 = (v12 * dpv2 - v1v2 * dpv1) * det;

		return (Magnitude(dp + v2 * t2 - v1 * t1));
	}

	Vector3D a = Cross(dp, v1);
	return (sqrt(Dot(a, a) / v12));
}

struct Plane
{
	float x, y, z, w;

	Plane() = default;

	Plane(float nx, float ny, float nz, float d) {
		x = nx;
		y = ny;
		z = nz;
		w = d;
	}

	Plane(const Vector3D& n, float d) {
		x = n.x;
		y = n.y;
		z = n.z;
		w = d;
	}

	const Vector3D& GetNormal(void) const {
		return (reinterpret_cast<const Vector3D&>(x));
	}

};

float Dot(const Plane& f, const Vector3D& v) {
	return (f.x * v.x + f.y * v.y + f.z * v.z);
}

float Dot(const Plane& f, const Point3D& p) {
	return (f.x * p.x + f.y * p.y + f.z * p.z);
}

struct Line
{
	Vector3D direction;
	Vector3D moment;

	Line() = default;

	Line(float vx, float vy, float vz, float mx, float my, float mz) :
		direction(vx, vy, vz), moment(mx, my, mz) {
		}

	Line(const Vector3D& v, const Vector3D& m) {
		direction = v;
		moment = m;
	}
};
	
int
main()
{
	std::cout << "First program\n";
	return 0;
}

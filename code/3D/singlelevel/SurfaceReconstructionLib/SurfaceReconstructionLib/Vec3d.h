#pragma once
#include <cmath>

struct Vec3d {
public:
	Vec3d() {
		x = 0; y = 0; z = 0;
	}
	Vec3d(double value1, double value2, double value3) {
		x = value1; y = value2; z = value3;
	}
	double x, y, z;
	inline Vec3d operator-() const;  // ���ظ���
	inline Vec3d operator+(const Vec3d& vec) const;  // �������
	inline Vec3d operator-(const Vec3d& vec) const;  // �������
	inline Vec3d operator*(double c) const;  // �˱���
	inline Vec3d operator/(double c) const;  // ������
	inline double operator*(const Vec3d& vec) const;  // �ڻ�
	inline Vec3d Cross(const Vec3d& vec) const;  // ���
	inline void operator+=(const Vec3d& vec);  // ���ظ�ֵ�����+=
	inline void operator-=(const Vec3d& vec);  // ���ظ�ֵ�����-=
	inline void operator*=(double c);  // ���ظ�ֵ�����*=
	inline void operator/=(double c);  // ���ظ�ֵ�����/=
	//inline void operator=(const Vec3d& vec);  // ���ظ�ֵ�����=
	inline double Length() const;  // ��ģ
	inline double SquaredLength() const;  // ��ģƽ��
	inline void Normalize();  // ��λ��
	inline Vec3d GetNormalizedVec() const;  // ���ɵ�λ������
	inline void ApplyMultiplication(const Vec3d& vec);  // ��Ԫ�����
	inline void ApplySqrt();  // ��ÿ��Ԫ�ؿ���
	inline void ApplyAbs();  // ��ÿ��Ԫ��ȡ����ֵ
};

// ���ظ���
inline Vec3d Vec3d::operator-() const { return Vec3d(-x, -y, -z); }
// �����Ӽ�
inline Vec3d Vec3d::operator+(const Vec3d& vec) const { return Vec3d(x + vec.x, y + vec.y, z + vec.z); }
inline Vec3d Vec3d::operator-(const Vec3d& vec) const { return Vec3d(x - vec.x, y - vec.y, z - vec.z); }
// ������ĳ˳�
inline Vec3d Vec3d::operator*(double c) const { return Vec3d(x*c, y*c, z*c); }
inline Vec3d Vec3d::operator/(double c) const { return Vec3d(x / c, y / c, z / c); }
// �ڻ������
inline double Vec3d::operator*(const Vec3d& vec) const { return x * vec.x + y * vec.y + z * vec.z; }
inline Vec3d Vec3d::Cross(const Vec3d& vec) const {
	return Vec3d(
		y*vec.z - vec.y*z,
		vec.x*z - x * vec.z,
		x*vec.y - vec.x*y
	);
}
// ���ظ�ֵ�����
inline void Vec3d::operator+=(const Vec3d& vec) { x += vec.x; y += vec.y; z += vec.z; }
inline void Vec3d::operator-=(const Vec3d& vec) { x -= vec.x; y -= vec.y; z -= vec.z; }
inline void Vec3d::operator*=(double c) { x *= c; y *= c; z *= c; }
inline void Vec3d::operator/=(double c) { x /= c; y /= c; z /= c; }
//inline void Vec3d::operator=(const Vec3d& vec) { x = vec.x; y = vec.y; z = vec.z; }
// �ṩ��ģ����ģƽ������λ������
inline double Vec3d::Length() const { return sqrt(x*x + y * y + z * z); }
inline double Vec3d::SquaredLength() const { return x * x + y * y + z * z; }
inline void Vec3d::Normalize() {
	double len = sqrt(x*x + y * y + z * z);
	x /= len; y /= len; z /= len;
}
inline Vec3d Vec3d::GetNormalizedVec() const {
	double len = sqrt(x*x + y * y + z * z);
	return Vec3d(x / len, y / len, z / len);
}
// Ϊÿһ��ֵӦ�ñ任
inline void Vec3d::ApplyMultiplication(const Vec3d& vec) {
	x *= vec.x; y *= vec.y; z *= vec.z;
}
inline void Vec3d::ApplySqrt() {
	x = sqrt(x); y = sqrt(y); z = sqrt(z);
};
inline void Vec3d::ApplyAbs() {
	x = abs(x); y = abs(y); z = abs(z);
}
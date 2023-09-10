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
	inline Vec3d operator-() const;  // 重载负号
	inline Vec3d operator+(const Vec3d& vec) const;  // 向量相加
	inline Vec3d operator-(const Vec3d& vec) const;  // 向量相减
	inline Vec3d operator*(double c) const;  // 乘标量
	inline Vec3d operator/(double c) const;  // 除标量
	inline double operator*(const Vec3d& vec) const;  // 内积
	inline Vec3d Cross(const Vec3d& vec) const;  // 外积
	inline void operator+=(const Vec3d& vec);  // 重载赋值运算符+=
	inline void operator-=(const Vec3d& vec);  // 重载赋值运算符-=
	inline void operator*=(double c);  // 重载赋值运算符*=
	inline void operator/=(double c);  // 重载赋值运算符/=
	//inline void operator=(const Vec3d& vec);  // 重载赋值运算符=
	inline double Length() const;  // 求模
	inline double SquaredLength() const;  // 求模平方
	inline void Normalize();  // 单位化
	inline Vec3d GetNormalizedVec() const;  // 生成单位化向量
	inline void ApplyMultiplication(const Vec3d& vec);  // 按元素相乘
	inline void ApplySqrt();  // 对每个元素开根
	inline void ApplyAbs();  // 对每个元素取绝对值
};

// 重载负号
inline Vec3d Vec3d::operator-() const { return Vec3d(-x, -y, -z); }
// 向量加减
inline Vec3d Vec3d::operator+(const Vec3d& vec) const { return Vec3d(x + vec.x, y + vec.y, z + vec.z); }
inline Vec3d Vec3d::operator-(const Vec3d& vec) const { return Vec3d(x - vec.x, y - vec.y, z - vec.z); }
// 与标量的乘除
inline Vec3d Vec3d::operator*(double c) const { return Vec3d(x*c, y*c, z*c); }
inline Vec3d Vec3d::operator/(double c) const { return Vec3d(x / c, y / c, z / c); }
// 内积和外积
inline double Vec3d::operator*(const Vec3d& vec) const { return x * vec.x + y * vec.y + z * vec.z; }
inline Vec3d Vec3d::Cross(const Vec3d& vec) const {
	return Vec3d(
		y*vec.z - vec.y*z,
		vec.x*z - x * vec.z,
		x*vec.y - vec.x*y
	);
}
// 重载赋值运算符
inline void Vec3d::operator+=(const Vec3d& vec) { x += vec.x; y += vec.y; z += vec.z; }
inline void Vec3d::operator-=(const Vec3d& vec) { x -= vec.x; y -= vec.y; z -= vec.z; }
inline void Vec3d::operator*=(double c) { x *= c; y *= c; z *= c; }
inline void Vec3d::operator/=(double c) { x /= c; y /= c; z /= c; }
//inline void Vec3d::operator=(const Vec3d& vec) { x = vec.x; y = vec.y; z = vec.z; }
// 提供求模、求模平方、单位化操作
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
// 为每一个值应用变换
inline void Vec3d::ApplyMultiplication(const Vec3d& vec) {
	x *= vec.x; y *= vec.y; z *= vec.z;
}
inline void Vec3d::ApplySqrt() {
	x = sqrt(x); y = sqrt(y); z = sqrt(z);
};
inline void Vec3d::ApplyAbs() {
	x = abs(x); y = abs(y); z = abs(z);
}
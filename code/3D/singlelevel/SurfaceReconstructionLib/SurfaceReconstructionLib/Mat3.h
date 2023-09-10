#pragma once
#include "Vec3d.h"
#define SWAP(T,a,b) {T t; t=a, a=b, b=t;}

struct Mat3 {
public:
	Mat3() {
		u = Vec3d();
		v = Vec3d();
		w = Vec3d();
	}
	Mat3(Vec3d u, Vec3d v, Vec3d w) {
		this->u = u;
		this->v = v;
		this->w = w;
	}
	Vec3d u, v, w;
	// 将矩阵列看做一组基，进行基变换
	inline Vec3d BaseTransform(Vec3d original) {
		return Vec3d(u*original, v*original, w*original);
	}
	// 使用列主元高斯消去法求解Ax=b
	Vec3d SolveLinerEquation(Vec3d b) {
		Mat3 A = Mat3(u, v, w);
		// 第一列主元交换与消元
		if ((abs(A.u.x) < abs(A.u.y)) || (abs(A.u.x) < abs(A.u.z))) {
			if (abs(A.u.z) < abs(A.u.y)) {
				SWAP(double, A.u.x, A.u.y);
				SWAP(double, A.v.x, A.v.y);
				SWAP(double, A.w.x, A.w.y);
				SWAP(double, b.x, b.y);
			}
			else {
				SWAP(double, A.u.x, A.u.z);
				SWAP(double, A.v.x, A.v.z);
				SWAP(double, A.w.x, A.w.z);
				SWAP(double, b.x, b.z);
			}
		}
		if (A.u.x > 1e-30) {  // 可解决欠定方程组问题
			A.v.y -= A.v.x * A.u.y / A.u.x;
			A.w.y -= A.w.x * A.u.y / A.u.x;
			b.y -= b.x * A.u.y / A.u.x;
			A.v.z -= A.v.x * A.u.z / A.u.x;
			A.w.z -= A.w.x * A.u.z / A.u.x;
			b.z -= b.x * A.u.z / A.u.x;
		}
		// 第二列主元交换与消元
		if (abs(A.v.y) < abs(A.v.z)) {
			SWAP(double, A.v.y, A.v.z);
			SWAP(double, A.w.y, A.w.z);
			SWAP(double, b.y, b.z);
		}
		if (A.v.y > 1e-30) {  // 可解决欠定方程组问题
			A.w.z -= A.w.y * A.v.z / A.v.y;
			b.z -= b.y * A.v.z / A.v.y;
		}
		// 回代求解
		Vec3d res;
		if (A.w.z > 1e-30) {
			res.z = b.z / A.w.z;
		} else { res.z = 0; }
		if (A.v.y > 1e-30) {
			res.y = (b.y - A.w.y*res.z) / A.v.y;
		} else { res.y = 0; }
		if (A.u.x > 1e-30) {
			res.x = (b.x - A.w.x*res.z - A.v.x*res.y) / A.u.x;
		} else { res.x = 0; }
		// 若系数矩阵秩为2，则解空间必定与至少一坐标平面相交，若系数矩阵秩为1，则则解空间必定与至少一坐标轴相交
		return res;
	}
};
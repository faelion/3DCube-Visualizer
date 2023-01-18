#pragma once

#include "Matrix.h"
#include "Quaternion.h"
#include "Vector4.h"
#include "Vector3.h"
#include "Vector2.h"
#include <math.h>

#define PI 3.14159265359
#define TORAD (PI / 180.0)
#define DEGTORAD(deg) deg * TORAD

inline Vector3 operator*(const Matrix<double, 3, 3>& mat, const Vector3& v) {
	Vector3 result = { 0, 0, 0 };

	result.x = v.x * mat[0][0] + v.y * mat[0][1] + v.z * mat[0][2];
	result.y = v.x * mat[1][0] + v.y * mat[1][1] + v.z * mat[1][2];
	result.z = v.x * mat[2][0] + v.y * mat[2][1] + v.z * mat[2][2];

	return result;
}

inline Vector4 operator*(const Matrix<double, 4, 4>& mat, const Vector4& v) {
	Vector4 result;

	result.x = v.x * mat[0][0] + v.y * mat[0][1] + v.z * mat[0][2];
	result.y = v.x * mat[1][0] + v.y * mat[1][1] + v.z * mat[1][2];
	result.z = v.x * mat[2][0] + v.y * mat[2][1] + v.z * mat[2][2];
	result.w = v.w * mat[3][0] + v.y * mat[3][1] + v.z * mat[3][2];

	return result;
}

namespace Math {
	enum class RotationType {
		ROTATION_MATRIX,
		EULER_ANGLES,
		EULER_AXIS_ANGLE,
		QUATERNION,
		ROTATION_VECTOR
	};

	struct EulerAxisAngle {
		Vector3 eulerAxis;
		double angle;
	};

	struct Rotations {
		Matrix<double, 3, 3> rotationMatrix;
		Vector3 eulerAngles;
		EulerAxisAngle eulerAxisAngle;
		Quaternion quaternion;
		Vector3 rotationVector;
	};

	inline Matrix<double, 2, 2> RotMatrix(double yaw) {
		Matrix<double, 2, 2> mat;

		mat[0][0] = cos(yaw);
		mat[0][1] = sin(yaw);
		mat[1][0] = -sin(yaw);
		mat[1][1] = cos(yaw);

		return mat;
	}

	inline Vector2 getTranslation(const Vector2& v1, const Vector2& v2, const Matrix<double, 2, 2>& rotMatrix) {
		Vector2 vec;

		vec.x = v1.x - rotMatrix[0][0] * v2.x - rotMatrix[0][1] * v2.y;
		vec.y = v1.y - rotMatrix[1][0] * v2.x - rotMatrix[1][1] * v2.y;

		return vec;
	}

	inline Vector3 Translate(const Vector3& vec, const Matrix<double, 4, 4>& mat) {
		// Homogeneous 3d vector
		Vector4 vecAug = vec.Homogeneous();

		// Dot product matrix * vector
		Vector4 resultAug = mat * vecAug;

		// Get back to vector2 (z = 1.0)
		Vector3 result = { resultAug.x, resultAug.y, resultAug.z };

		return result;
	}

	inline Matrix<double, 3, 3> AugmentMatrix(const Matrix<double, 2, 2>& mat, const Vector2& vec) {
		Matrix<double, 3, 3> matrix;

		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				matrix[i][j] = mat[i][j];
			}
		}

		matrix[0][2] = vec.x;
		matrix[1][2] = vec.y;
		matrix[2][2] = 1.0;

		return matrix;
	}

	inline Matrix<double, 4, 4> AugmentMatrix(const Matrix<double, 3, 3>& mat, const Vector3& vec) {
		Matrix<double, 4, 4> matrix;

		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				matrix[i][j] = mat[i][j];
			}
		}

		matrix[0][3] = vec.x;
		matrix[1][3] = vec.y;
		matrix[2][3] = vec.z;
		matrix[3][3] = 1.0;

		return matrix;
	}

	template<class T, int R, int C>
	inline Matrix<T, R, C> EulerToRotMatrix(double roll, double pitch, double yaw) {
		Matrix<T, R, C> rotMatrix;

		double cp = cos(pitch);
		double sp = sin(pitch);

		double cy = cos(yaw);
		double sy = sin(yaw);
		
		double cr = cos(roll);
		double sr = sin(roll);

		rotMatrix[0][0] = cy * cp;
		rotMatrix[0][1] = cy * sp * sr - sy * cr;
		rotMatrix[0][2] = cy * sp * cr + sy * sr;

		rotMatrix[1][0] = sy * cp;
		rotMatrix[1][1] = sy * sp * sr + cy * cr;
		rotMatrix[1][2] = sy * sp * cr - cy * sr;

		rotMatrix[2][0] = -sp;
		rotMatrix[2][1] = cp * sr;
		rotMatrix[2][2] = cp * cr;

		return rotMatrix;
	}

	template<class T, int R, int C>
	inline Matrix<T, R, C> EulerToRotMatrix(Vector3 eulers) {
		return EulerToRotMatrix<T, R, C>(eulers.x, eulers.y, eulers.z);
	}

	template<class T, int R, int C>
	inline Matrix<T, R, C> EulerAxisAngleToRotMatrix(Vector3 axis, double angle) {
		Matrix<T, R, C> rotMatrix;

		axis.Normalize();

		double c = cos(angle);
		double s = sin(angle);
		double t = 1.0 - c;

		rotMatrix[0][0] = c + axis.x * axis.x * t;
		rotMatrix[1][1] = c + axis.y * axis.y * t;
		rotMatrix[2][2] = c + axis.z * axis.z * t;

		double tmp1 = axis.x * axis.y * t;
		double tmp2 = axis.z * s;

		rotMatrix[1][0] = tmp1 + tmp2;
		rotMatrix[0][1] = tmp1 - tmp2;

		tmp1 = axis.x * axis.z * t;
		tmp2 = axis.y * s;

		rotMatrix[2][0] = tmp1 - tmp2;
		rotMatrix[0][2] = tmp1 + tmp2;

		tmp1 = axis.y * axis.z * t;
		tmp2 = axis.x * s;

		rotMatrix[2][1] = tmp1 + tmp2;
		rotMatrix[1][2] = tmp1 - tmp2;

		return rotMatrix;
	}

	inline EulerAxisAngle EulerToEulerAxisAngle(Vector3 eulerAngles) {
		EulerAxisAngle eulerAxisAngle;

		double cx = cos(eulerAngles.x / 2);
		double sx = sin(eulerAngles.x / 2);

		double cy = cos(eulerAngles.y / 2);
		double sy = sin(eulerAngles.y / 2);

		double cz = cos(eulerAngles.z / 2);
		double sz = sin(eulerAngles.z / 2);

		double cxcy = cx * cy;
		double sxsy = sx * sy;

		eulerAngles.x = cxcy * sz + sxsy * cz;
		eulerAngles.y = sx * cy * cz + cx * sy * sz;
		eulerAngles.z = cx * sy * cz - sx * cy * sz;

		eulerAxisAngle.angle = 2 * acos(cxcy * cz - sxsy * sz);

		eulerAngles.Normalize();

		eulerAxisAngle.eulerAxis = eulerAngles;

		return eulerAxisAngle;
	}

	inline Quaternion EulerToQuaternion(Vector3 eulerAngles) {
		Quaternion quat;

		double cx = cos(eulerAngles.x / 2);
		double sx = sin(eulerAngles.x / 2);

		double cy = cos(eulerAngles.y / 2);
		double sy = sin(eulerAngles.y / 2);

		double cz = cos(eulerAngles.z / 2);
		double sz = sin(eulerAngles.z / 2);

		double cxcy = cx * cy;
		double sxsy = sx * sy;

		quat.x = cxcy * sz + sxsy * cz;
		quat.y = sx * cy * cz + cx * sy * sz;
		quat.z = cx * sy * cz - sx * cy * sz;

		quat.w = cxcy * cz - sxsy * sz;

		return quat;
	}

	template<class T, int R, int C>
	inline Vector3 RotMatrixToEuler(const Matrix<T, R, C>& mat) {
		Vector3 v;

		v.x = atan2(mat[2][1], mat[2][2]);
		v.y = atan2(-mat[2][0], sqrt(mat[2][1] * mat[2][1] + mat[2][2] * mat[2][2]));
		v.z = atan2(mat[1][0], mat[0][0]);

		return v;
	}

	template<class T, int R, int C>
	inline EulerAxisAngle RotMatrixToEulerAxisAngles(const Matrix<T, R, C>& mat) {
		EulerAxisAngle eulerAxisAngle;

		Matrix<T, R, C> transposed = mat.Transpose();
		Matrix<T, R, C> mvp;

		eulerAxisAngle.angle = acos((mat.Trace() - 1) / 2);

		double as = 2 * sin(eulerAxisAngle.angle);

		if (!as) as = 1;

		for (int i = 0; i < R; i++)
		{
			for (int j = 0; j < C; j++)
			{
				mvp[i][j] = (mat[i][j] - transposed[i][j]) / as;
			}
		}

		eulerAxisAngle.eulerAxis.x = (mvp[2][1] - mvp[1][2]) / 2;
		eulerAxisAngle.eulerAxis.y = (mvp[0][2] - mvp[2][0]) / 2;
		eulerAxisAngle.eulerAxis.z = (mvp[1][0] - mvp[0][1]) / 2;

		return eulerAxisAngle;
	}

	template<class T, int R, int C>
	inline Quaternion RotMatrixToQuaternion(const Matrix<T, R, C>& mat) {
		Quaternion quat;

		/*double trace = mat[0][0] + mat[1][1] + mat[2][2];
		if (trace > 0) {
			float s = 0.5f / sqrt(trace + 1.0f);
			quat.w = 0.25f / s;
			quat.x = (mat[2][1] - mat[1][2]) * s;
			quat.y = (mat[0][2] - mat[2][0]) * s;
			quat.z = (mat[1][0] - mat[0][1]) * s;
		}
		else {
			if (mat[0][0] > mat[1][1] && mat[0][0] > mat[2][2]) {
				float s = 2.0f * sqrt(1.0f + mat[0][0] - mat[1][1] - mat[2][2]);
				quat.w = (mat[1][2] - mat[2][1]) / s;
				quat.x = 0.25f * s;
				quat.y = (mat[0][1] + mat[1][0]) / s;
				quat.z = (mat[0][2] + mat[2][0]) / s;
			}
			else if (mat[1][1] > mat[2][2]) {
				float s = 2.0f * sqrt(1.0f + mat[1][1] - mat[0][0] - mat[2][2]);
				quat.w = (mat[0][2] - mat[2][0]) / s;
				quat.x = (mat[0][1] + mat[1][0]) / s;
				quat.y = 0.25f * s;
				quat.z = (mat[1][2] + mat[2][1]) / s;
			}
			else {
				float s = 2.0f * sqrt(1.0f + mat[2][2] - mat[0][0] - mat[1][1]);
				quat.w = (mat[0][1] - mat[1][0]) / s;
				quat.x = (mat[0][2] + mat[2][0]) / s;
				quat.y = (mat[1][2] + mat[2][1]) / s;
				quat.z = 0.25f * s;
			}
		}*/

		quat.w = sqrt(1.0 + mat[0][0] + mat[1][1] + mat[2][2]) / 2;
		double w4 = quat.w * 4.0;

		if (!w4) w4 = 1;

		quat.x = (mat[2][1] - mat[1][2]) / w4;
		quat.y = (mat[0][2] - mat[2][0]) / w4;
		quat.z = (mat[1][0] - mat[0][1]) / w4;

		return quat;
	}

	inline Quaternion EulerAxisAngleToQuaternion(Vector3 axis, double angle)
	{
		Quaternion quat;

		double s = sin(angle / 2);

		quat.x = axis.x * s;
		quat.y = axis.y * s;
		quat.z = axis.z * s;
		quat.w = cos(angle / 2);

		return quat;
	}

	inline Quaternion EulerAxisAngleToQuaternion(EulerAxisAngle eulerAxisAngle) {
		return EulerAxisAngleToQuaternion(eulerAxisAngle.eulerAxis, eulerAxisAngle.angle);
	}

	inline Vector3 EulerAxisAngleToEulerAngles(EulerAxisAngle eulerAxisAngle) {
		Vector3 eulerAngles;
		Vector3 eulerAxis = eulerAxisAngle.eulerAxis;

		eulerAxis.Normalize();

		double c = cos(eulerAxisAngle.angle);
		double s = sin(eulerAxisAngle.angle);
		double t = 1.0 - c;

		if ((eulerAxis.x*eulerAxis.y*t + eulerAxis.z * s) > 0.998) { // north pole singularity detected
			eulerAngles.x = 2 * atan2(eulerAxis.x*sin(eulerAxisAngle.angle / 2), cos(eulerAxisAngle.angle / 2));
			eulerAngles.y = PI / 2;
			eulerAngles.z = 0;
		}else if ((eulerAxis.x*eulerAxis.y*t + eulerAxis.z * s) < -0.998) { // south pole singularity detected
			eulerAngles.x = -2 * atan2(eulerAxis.x*sin(eulerAxisAngle.angle / 2), cos(eulerAxisAngle.angle / 2));
			eulerAngles.y = -PI / 2;
			eulerAngles.z = 0;
		}
		else {
			eulerAngles.x = atan2(eulerAxis.y * s - eulerAxis.x * eulerAxis.z * t, 1 - (eulerAxis.y * eulerAxis.y + eulerAxis.z * eulerAxis.z) * t);
			eulerAngles.y = asin(eulerAxis.x * eulerAxis.y * t + eulerAxis.z * s);
			eulerAngles.z = atan2(eulerAxis.x * s - eulerAxis.y * eulerAxis.z * t, 1 - (eulerAxis.x * eulerAxis.x + eulerAxis.z * eulerAxis.z) * t);
		}

		return eulerAngles;
	}

	template<class T, int R, int C>
	inline Matrix<T, R, C> EulerAxisAngleToRotMatrix(EulerAxisAngle eulerAxisAngle) {
		Matrix<T, R, C> rotMatrix;

		Vector3 eulerAxis = eulerAxisAngle.eulerAxis;

		double c = cos(eulerAxisAngle.angle);
		double s = sin(eulerAxisAngle.angle);
		double t = 1.0 - c;

		eulerAxis.Normalize();

		rotMatrix[0][0] = c + eulerAxis.x * eulerAxis.x * t;
		rotMatrix[1][1] = c + eulerAxis.y * eulerAxis.y * t;
		rotMatrix[2][2] = c + eulerAxis.z * eulerAxis.z * t;

		double tmp1 = eulerAxis.x * eulerAxis.y * t;
		double tmp2 = eulerAxis.z * s;

		rotMatrix[1][0] = tmp1 + tmp2;
		rotMatrix[0][1] = tmp1 - tmp2;

		tmp1 = eulerAxis.x * eulerAxis.z * t;
		tmp2 = eulerAxis.y * s;

		rotMatrix[2][0] = tmp1 - tmp2;
		rotMatrix[0][2] = tmp1 + tmp2;

		tmp1 = eulerAxis.y * eulerAxis.z * t;
		tmp2 = eulerAxis.x * s;

		rotMatrix[2][1] = tmp1 + tmp2;
		rotMatrix[1][2] = tmp1 - tmp2;

		return rotMatrix;
	}

	inline Vector3 EulerAxisAngleToRotVector(EulerAxisAngle eulerAxisAngle) {
		return eulerAxisAngle.eulerAxis * eulerAxisAngle.angle;
	}

	inline Vector3 QuaternionToEuler(Quaternion quat) {
		Vector3 eulerAngles;

		quat.Normalize();

		double pole = quat.x * quat.y + quat.z * quat.w;

		if (pole > 0.499) { // North pole
			eulerAngles.x = 2 * atan2(quat.x, quat.w);
			eulerAngles.y = PI / 2;
			eulerAngles.z = 0;
		}
		else if (pole < -0.499) { // South pole
			eulerAngles.x = -2 * atan2(quat.x, quat.w);
			eulerAngles.y = -PI / 2;
			eulerAngles.z = 0;
		}
		else {
			double sqx = quat.x*quat.x;
			double sqy = quat.y*quat.y;
			double sqz = quat.z*quat.z;
			eulerAngles.x = atan2(2 * quat.y*quat.w - 2 * quat.x*quat.z, 1 - 2 * sqy - 2 * sqz);
			eulerAngles.y = asin(2 * pole);
			eulerAngles.z = atan2(2 * quat.x*quat.w - 2 * quat.y*quat.z, 1 - 2 * sqx - 2 * sqz);
		}

		return eulerAngles;
	}

	inline EulerAxisAngle QuaternionToEulerAxisAngle(Quaternion quat) { 
		EulerAxisAngle eulerAxisAngle;

		if (quat.w > 1) quat.Normalize();

		eulerAxisAngle.angle = 2 * acos(quat.w);

		double s = sqrt(1 - quat.w * quat.w);

		if (s != 0) {
			eulerAxisAngle.eulerAxis.x = quat.x / s;
			eulerAxisAngle.eulerAxis.y = quat.y / s;
			eulerAxisAngle.eulerAxis.z = quat.z / s;
		}
		else {
			eulerAxisAngle.eulerAxis.x = quat.x;
			eulerAxisAngle.eulerAxis.y = quat.y;
			eulerAxisAngle.eulerAxis.z = quat.z;
		}

		return eulerAxisAngle;
	}

	template<class T, int R, int C>
	inline Matrix<T, R, C> QuaternionToRotMatrix(Quaternion quat) {
		Matrix<T, R, C> rotMatrix;

		quat.Normalize();

		double sqx = quat.x*quat.x;
		double sqy = quat.y*quat.y;
		double sqz = quat.z*quat.z;
		double sqw = quat.w*quat.w;

		rotMatrix[0][0] = sqx - sqy - sqz + sqw;
		rotMatrix[1][1] = -sqx + sqy - sqz + sqw;
		rotMatrix[2][2] = -sqx - sqy + sqz + sqw;

		double tmp1 = quat.x * quat.y;
		double tmp2 = quat.z * quat.w;

		rotMatrix[1][0] = 2.0 * (tmp1 + tmp2);
		rotMatrix[0][1] = 2.0 * (tmp1 - tmp2);

		tmp1 = quat.x * quat.z;
		tmp2 = quat.y * quat.w;
		rotMatrix[2][0] = 2.0 * (tmp1 - tmp2);
		rotMatrix[0][2] = 2.0 * (tmp1 + tmp2);

		tmp1 = quat.y * quat.z;
		tmp2 = quat.x * quat.w;
		rotMatrix[2][1] = 2.0 * (tmp1 + tmp2);
		rotMatrix[1][2] = 2.0 * (tmp1 - tmp2);

		return rotMatrix;
	}

	inline EulerAxisAngle RotVectorToEulerAxisAngle(Vector3 rotVector) {
		EulerAxisAngle eulerAxisAngle;

		eulerAxisAngle.angle = sqrt(rotVector.x * rotVector.x + rotVector.y * rotVector.y + rotVector.z * rotVector.z);

		if (eulerAxisAngle.angle) eulerAxisAngle.eulerAxis = rotVector / eulerAxisAngle.angle;
		else eulerAxisAngle.eulerAxis = rotVector;

		return eulerAxisAngle;
	}

	inline void CalculateRotations(Rotations& rotations, RotationType rotationType) {
		switch (rotationType) {
			case RotationType::EULER_ANGLES:
				rotations.rotationMatrix = EulerToRotMatrix<double, 3, 3>(rotations.eulerAngles);
				break;
			case RotationType::EULER_AXIS_ANGLE:
				rotations.rotationMatrix = EulerAxisAngleToRotMatrix<double, 3, 3>(rotations.eulerAxisAngle);
				break;
			case RotationType::QUATERNION:
				rotations.rotationMatrix = QuaternionToRotMatrix<double, 3, 3>(rotations.quaternion);
				break;
			case RotationType::ROTATION_MATRIX:
				break;
			case RotationType::ROTATION_VECTOR:
				rotations.eulerAxisAngle = RotVectorToEulerAxisAngle(rotations.rotationVector);
				rotations.rotationMatrix = EulerAxisAngleToRotMatrix<double, 3, 3>(rotations.eulerAxisAngle);
				break;
		}

		rotations.eulerAngles = RotMatrixToEuler(rotations.rotationMatrix);
		rotations.eulerAxisAngle = RotMatrixToEulerAxisAngles(rotations.rotationMatrix);
		rotations.quaternion = RotMatrixToQuaternion(rotations.rotationMatrix);
		rotations.rotationVector = EulerAxisAngleToRotVector(rotations.eulerAxisAngle);
	}
}
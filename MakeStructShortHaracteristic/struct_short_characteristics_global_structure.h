#pragma once
#ifndef SHORT_CHARACTERISTICS_GLOBAL_H
#define SHORT_CHARACTERISTICS_GLOBAL_H

#include "struct_short_characteristics_headers.h"

#define BASE_ADRESS "D:\\Desktop\\FilesCourse\\IllumGrid\\"


#ifdef _MSC_VER
#define fwrite_unlocked _fwrite_nolock
#define fread_unlocked  _fread_nolock
#endif

typedef int IntId;
typedef double Type;
typedef std::string Str_Type;

typedef Eigen::Vector3d Vector3;
typedef Eigen::Vector2d Vector2;
typedef Eigen::Matrix3d Matrix3;

extern std::vector<Vector3> X;
extern std::vector<Vector2> X0;
extern std::vector<Type> S;

using namespace std;
using namespace std::chrono;

struct Normals {
	std::vector<Vector3> n;
	Normals() {
	}

	Normals(const int size) {
		n.resize(size);
	}
};

struct cell {
	//int id;  - ����� � �������
	std::vector<Vector3> nodes_value;
	std::vector<int> neighbours_id_face;

	cell() {
		//id = -1;
		nodes_value.resize(4, Vector3(-666, -666, -666));
		neighbours_id_face.resize(4, -1);
	}
};

#define PI 3.14159265358979323846
const double eps = 1e-10;

extern Vector3 start_point_plane_coord;   // ������ ��������� ���������
extern Matrix3 transform_matrix;          // ������� �������� �� �������� ��������� � ���������
extern Matrix3 inverse_transform_matrix;  // ������� �������� �� ��������� � ������� ��������

extern Matrix3	straight_face;  // 3 ���� ������������
extern Matrix3 inclined_face;  // 3 ���� ������������ �� ��������� ���������

//std::vector<Type> Illum2;

extern const Vector3 center_point;
extern const Type inner_radius;

// ��������� ������ ����� (unstructured_grid)
extern vtkDataArray* density;
extern vtkDataArray* absorp_coef;
extern vtkDataArray* rad_en_loose_rate;

extern Type square_surface;  // ������� ����������� ���������� 

extern Vector3 center_local_sphere;  // ����� ��������� ����� ����� ������������ ���������

extern int num_cur_direction; // ����� �������� �����������
extern Vector3 cur_direction;

extern int count_negative_interpolation; // ����� ������������� �������� ������������ ������������


// ��������� ����� � ���������� �����:
const Type Rsphere = 0.1;
const Type R1disk = 0.1;
const Type R2disk = 0.2;

extern std::vector<Type> res_inner_bound;  // �������� �� ���������� �������
extern std::vector<Vector3> x_try_surface;
extern std::vector<int> id_try_surface;

extern int pos_x_try;
extern int posX;
extern int posX0;
extern int posOutC;
extern int posOut;
extern int posIn;
extern int posS;
extern int posRes;

extern int n_out;


#endif
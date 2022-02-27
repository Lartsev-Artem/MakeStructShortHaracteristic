#include "struct_short_characteristics_main.h"
Vector3 start_point_plane_coord;   // начало координат плоскости
Matrix3 transform_matrix;          // матрица перехода из базового тетраэдра в плоскость
Matrix3 inverse_transform_matrix;  // матрица перехода из плоскости в базовый тетраэдр

Matrix3	straight_face;  // 3 узла интерпол€ции
Matrix3 inclined_face;  // 3 узла интерпол€ции на наклонной плоскости

//std::vector<Type> Illum2;

// скал€рные данные сетки (unstructured_grid)
vtkDataArray* density;
vtkDataArray* absorp_coef;
vtkDataArray* rad_en_loose_rate;

Type square_surface;  // площадь поверхности дискретной 

Vector3 center_local_sphere;  // центр описанной сферы около стандартного тетраэдра

int num_cur_direction; // номер текущего направлени€
Vector3 cur_direction;

int count_negative_interpolation; // число отрицательных значений интерпол€ции интесивности

std::vector<Vector3> X;
std::vector<Vector2> X0;
std::vector<Type> S;

std::vector<Type> res_inner_bound;
std::vector<Vector3> x_try_surface;
std::vector<int> id_try_surface;
int pos_x_try;

std::vector<int> in_id;
std::vector<int> out_id;

const Vector3 center_point(0, 0, 0);
const Type inner_radius = 0.51; // радиус внутренней сферы (с запасом)

//const Vector3 center_point(1, 0, 0);
//const Type R = 0.11;  // радиус внутренней сферы (с запасом)

int main(int argc, char* argv[])
{
	std::string name_file_settings = "";

	if (argc <= 1)
		name_file_settings = "D:\\Desktop\\FilesCourse\\settings_file_make.txt";
	else
		name_file_settings = argv[1];


	std::string name_file_vtk;
	std::string name_file_sphere_direction;
	std::string name_file_graph;
	std::string out_file_grid_vtk;
	std::string name_file_normals;

	if (ReadStartSettings(name_file_settings, name_file_vtk, name_file_sphere_direction, out_file_grid_vtk, name_file_graph, name_file_normals)) {
		std::cout << "Error reading the start settings\n";
		return 1;
	}

	vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid =
		vtkSmartPointer<vtkUnstructuredGrid>::New();


	Type _clock = -omp_get_wtime();
	if (ReadFileVtk(0, name_file_vtk, unstructured_grid, density, absorp_coef, rad_en_loose_rate, true)) {
		std::cout << "Error reading the file vtk\n";
		return 1;
	}
	_clock += omp_get_wtime();
	std::cout << "\n Reading time of the vtk_grid file: " << _clock << "\n";


	vector<Vector3> directions;
	vector<Type> squares;

	_clock = -omp_get_wtime();
	ReadSphereDirectionDecartToSpherical(name_file_sphere_direction, directions, squares, square_surface);
	_clock += omp_get_wtime();
	std::cout << "\n Reading time of the sphere_direction file: " << _clock << "\n";

	std::vector<cell> nodes_value;
	{
		std::vector<int> all_pairs_face;
		/*_clock = -omp_get_wtime();
		int count_unique_face = FindNeighborsPairFace(unstructured_grid, all_pairs_face);

		InitNodesValue(all_pairs_face, nodes_value);

		_clock += omp_get_wtime();*/

		FindNeighborsPairFaceAndBoundaries(unstructured_grid, all_pairs_face);

		InitNodesValue(all_pairs_face, nodes_value);

		std::unique_ptr<FILE, int(*)(FILE*)> file_id_neighbors(fopen((std::string(BASE_ADRESS) + "id_neighbors.bin").c_str(), "wb"), fclose);
		if (!file_id_neighbors) printf("id_neighbors not open\n");

		fwrite_unlocked(all_pairs_face.data(), sizeof(int), all_pairs_face.size(), file_id_neighbors.get());

		fclose(file_id_neighbors.get());
	}
	std::cout << "\n Finding time of the all_pairs_face: " << _clock << "\n";


	const int count_directions = directions.size();
	const int count_cells = unstructured_grid->GetNumberOfCells();

	InitGlobalValue(start_point_plane_coord, transform_matrix, inverse_transform_matrix, straight_face, inclined_face);

	// ”пор€доченные индексы €чеек по данному направлению
	vector<IntId> sorted_id_cell(count_cells);

	Eigen::Matrix4d vertex_tetra;
	/* x1 x2 x3 x4
	*  y1 y2 y3 y4
	*  z1 z2 z3 z4
	*  1  1  1  1
	*/

	std::vector<Type> Illum;// (count_directions * count_cells, 0);
	std::vector<Type> Illum2;// (count_directions * count_cells, 0);


	std::vector<Normals> normals;
	ReadNormalFile(name_file_normals, normals);

	{
		int n;
		std::unique_ptr<FILE, int(*)(FILE*)> file_id(fopen((std::string(BASE_ADRESS) + "id_defining_faces.bin").c_str(), "rb"), fclose);
		if (!file_id) { printf("Error file_id\n"); return 1; }

		fread_unlocked(&n, sizeof(int), 1, file_id.get());
		id_try_surface.resize(n);
		fread_unlocked(id_try_surface.data(), sizeof(int), n, file_id.get());

		fclose(file_id.get());
	}

	//std::vector<Vector3> centers;
	//FindAllCenterOfTetra(unstructured_grid, centers);
	//{
	//	ofile_centers.open(std::string(BASE_ADRESS) + "centers.txt");

	//	if (!ofile_centers.is_open()) printf("centers not open\n");
	//	
	//	ofile_centers << centers.size() << '\n';		
	//	for (auto &el : centers)
	//	{
	//		ofile_centers << el[0] << el[1] << el[2] << '\n';
	//	}
	//	ofile_centers.close();
	//}

	int count = 0;
	Type norm = 0;
	Vector3 direction;

	ofstream ofile;
	ofile.open("File_with_Logs.txt");

	std::vector<int> nCount(count_directions * count_cells, 0);

	 X.reserve(3 * count_cells);
	 X0.reserve(3 * count_cells);
	 S.reserve(3 * count_cells);
	 res_inner_bound.reserve(1000); // размер внутренней границе
	
	_clock = -omp_get_wtime();
		{

			/*---------------------------------- далее FOR по направлени€м----------------------------------*/
			for (int num_direction = 0; num_direction < count_directions; ++num_direction)
			{				
				direction = directions[num_direction];

				num_cur_direction = num_direction;
				cur_direction = direction;

				ReadGraphBin(name_file_graph + to_string(num_direction) + ".bin", sorted_id_cell);
				//ResetNodesValue(nodes_value);

				Vector3 x;
				Vector3 x0;
				int num_cell;
				const int count_cells = unstructured_grid->GetNumberOfCells();
				X.clear();
				X0.clear();
				S.clear();
				res_inner_bound.clear();

				X.reserve(3 * count_cells);
				X0.reserve(3 * count_cells);
				S.reserve(3 * count_cells);
				res_inner_bound.reserve(1000); // размер внутренней границе
				

				std::unique_ptr<FILE, int(*)(FILE*)>file_out_id(fopen((std::string(BASE_ADRESS) + "OutId" + to_string(num_direction) + ".bin").c_str(), "wb"), fclose);
				std::unique_ptr<FILE, int(*)(FILE*)>file_in_id (fopen((std::string(BASE_ADRESS) + "InId" + to_string(num_direction) + ".bin").c_str(), "wb"), fclose);
				std::unique_ptr<FILE, int(*)(FILE*)>file_x(fopen((std::string(BASE_ADRESS) + "X" + to_string(num_direction) + ".bin").c_str(), "wb"), fclose);
				std::unique_ptr<FILE, int(*)(FILE*)>file_x0_local(fopen((std::string(BASE_ADRESS) + "locX0" + to_string(num_direction) + ".bin").c_str(), "wb"), fclose);
				std::unique_ptr<FILE, int(*)(FILE*)>file_s(fopen((std::string(BASE_ADRESS) + "S" + to_string(num_direction) + ".bin").c_str(), "wb"), fclose);
				std::unique_ptr<FILE, int(*)(FILE*)>file_res_bound(fopen((std::string(BASE_ADRESS) + "ResBound" + to_string(num_direction) + ".bin").c_str(), "wb"), fclose);
				
				std::unique_ptr<FILE, int(*)(FILE*)>file_x_defining_faces(fopen((std::string(BASE_ADRESS) + "x_defining_faces" + to_string(num_direction) + ".bin").c_str(), "rb"), fclose);
				if (!file_x_defining_faces) printf("file_x_defining_faces not open\n");

				int n = -1;				
				fread_unlocked(&n, sizeof(int), 1, file_x_defining_faces.get());
				x_try_surface.resize(n);
				fread_unlocked(x_try_surface.data(), sizeof(Vector3), n, file_x_defining_faces.get());
				fclose(file_x_defining_faces.get());
				pos_x_try = 0;



			    if (!file_out_id) printf("Outid not open\n");
				if (!file_in_id) printf("Inid not open\n");
				if (!file_x) printf("Xid not open\n");
				if (!file_x0_local) printf("X0id not open\n");
				if (!file_s) printf("Sid not open\n");
				
				

				int face_state[4];  // 0=> выход€ща€ грань,  1=> вход€ща€   

				/*---------------------------------- далее FOR по €чейкам----------------------------------*/
				for (int h = 0; h < count_cells; ++h) {
					num_cell = sorted_id_cell[h];

					SetVertexMatrix(num_cell, unstructured_grid, vertex_tetra);
					FindInAndOutFaces(direction, num_cell, normals, face_state);

					in_id.clear();
					out_id.clear();


					for (size_t num_out_face = 0; num_out_face < 4; ++num_out_face) {
						if (!face_state[num_out_face]) {  // выход€щие грани
							GetNodes(num_cell, unstructured_grid, unstructured_grid->GetCell(num_cell), num_out_face, vertex_tetra, face_state, direction,
								nodes_value, Illum2, directions, squares);
							out_id.push_back(num_out_face);
							//ofile_out_id << num_out_face << ' ';
						}
					}

					int n_out = out_id.size();
					int n_in = in_id.size();

					//ofile_count_out_id << n_out << '\n';  
					
					// ввести массив дл€ идентичного хранени€
					{
						nCount[count_directions * num_cell + num_direction] = n_out;  // мимо кеша
					}

					fwrite_unlocked(out_id.data(), sizeof(int), n_out, file_out_id.get());

					fwrite_unlocked(in_id.data(), sizeof(int), n_in, file_in_id.get());
				

					//Type I_k_dir = GetValueInCenterCell(num_cell, unstructured_grid, unstructured_grid->GetCell(num_cell), centers[num_cell], direction, vertex_tetra,
						//nodes_value, Illum2, directions, squares);

					//Illum[num_direction * count_cells + num_cell] = 0;// I_k_dir;
				}
				/*---------------------------------- конец FOR по €чейкам----------------------------------*/


				for (size_t i = 0; i < X.size(); i++)
				{
					fwrite_unlocked(X[i].data(), sizeof(Type), 3, file_x.get());
				}
				for (size_t i = 0; i < X0.size(); i++)
				{
					fwrite_unlocked(X0[i].data(), sizeof(Type), 2, file_x0_local.get());
				}

				fwrite_unlocked(S.data(), sizeof(Type), S.size(), file_s.get());

				int size = res_inner_bound.size();
				fwrite_unlocked(&size, sizeof(int), 1, file_res_bound.get());
				fwrite_unlocked(res_inner_bound.data(), sizeof(Type), size, file_res_bound.get());
				
				fclose(file_out_id.get());
				fclose(file_in_id.get());

				fclose(file_x.get());
				fclose(file_x0_local.get());
				fclose(file_s.get());

				fclose(file_res_bound.get());

								
				printf("End direction number #%d\n", num_direction);
			}
			/*---------------------------------- конец FOR по направлени€м----------------------------------*/
			
		}

		std::unique_ptr<FILE, int(*)(FILE*)>file_count_out_id(fopen((std::string(BASE_ADRESS) + "CountOutId"  + ".bin").c_str(), "wb"), fclose);
		fwrite_unlocked(nCount.data(), sizeof(int), nCount.size(), file_count_out_id.get());
		fclose(file_count_out_id.get());

		_clock += omp_get_wtime();
		//count++;
		//norm = NormIllum(Illum, Illum2);
		//std::cout << "Error:= " << norm << '\n';
		std::cout << "Time of iter: " << _clock << '\n';
	
		
		ofile << "Time of iter: " << _clock << '\n';		



	//Illum.swap(Illum2);
	ofile.close();

	/*vector<Type> energy(count_cells);
	MakeEnergy(Illum, squares, square_surface, energy);

	WriteFileSolution(out_file_grid_vtk, Illum, energy, unstructured_grid);*/
	return 0;
}



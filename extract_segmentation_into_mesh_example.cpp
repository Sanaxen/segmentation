#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/mesh_segmentation.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <random>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> SM;
typedef boost::graph_traits<SM>::face_descriptor face_descriptor;

#define VEC3_ADD(d,v,w)		(d)[0]=(v)[0]+(w)[0]; (d)[1]=(v)[1]+(w)[1]; (d)[2]=(v)[2]+(w)[2]
#define VEC3_SUB(d,v,w)		(d)[0]=(v)[0]-(w)[0]; (d)[1]=(v)[1]-(w)[1]; (d)[2]=(v)[2]-(w)[2]
#define VEC3_CROSS(d,v,w)	(d)[0]=(v)[1]*(w)[2] - (v)[2]*(w)[1]; \
							(d)[1]=(v)[2]*(w)[0] - (v)[0]*(w)[2]; \
							(d)[2]=(v)[0]*(w)[1] - (v)[1]*(w)[0]
#define VEC3_NORMALIZE(v)	{ float n=sqrt((v)[0]*(v)[0]+(v)[1]*(v)[1]+(v)[2]*(v)[2]); \
							if(fabs(n)>1e-6) { float m=1.0f/n; (v)[0]*=m; (v)[1]*=m; (v)[2]*=m; } }

void load_obj_to_off(const char *filename)
{
	unsigned int uiMesh;
	FILE *pFile = fopen(filename, "r");
	double *pVertices, *pNormals;
	unsigned int *pIndices;
	int i, j, f, iNumVertices, iNumIndices;
	unsigned int i0, i1, i2;
	double v1[3], v2[3], n[3];
	char buf[256];

	char tmpfile[256];
	FILE* fp;
	char drive_[_MAX_DRIVE];	// ドライブ名
	char dir_[_MAX_DIR];		// ディレクトリ名
	char fname_[_MAX_FNAME];	// ファイル名
	char ext_[_MAX_EXT];		// 拡張子

	if (!pFile)
	{
		fprintf(stderr, "cannot open file!\n");
		return;
	}

	_splitpath(filename, drive_, dir_, fname_, ext_);

	iNumVertices = 0;
	iNumIndices = 0;
	while (fgets(buf, 256, pFile) != NULL)
	{
		if (buf[0] == 'v' && !(buf[1] == 'n' || buf[1] == 't')) iNumVertices++;
		if (buf[0] == 'f') iNumIndices++;
	}
	fclose(pFile);

	pFile = fopen(filename, "r");

	printf("vertex %d face %d\n", iNumVertices, iNumIndices);
	/* read vertex data */
	iNumVertices *= 3;
	iNumIndices *= 3;
	pVertices = (double*)malloc(iNumVertices * sizeof(double));
	pNormals = (double*)calloc(iNumVertices, sizeof(double));
	pIndices = (unsigned int*)malloc(iNumIndices * sizeof(unsigned int));

	for (i = 0; i < iNumVertices / 3; ++i)
	{
		double x, y, z;
		fgets(buf, 256, pFile);
		while (buf[0] != 'v' || (buf[0] == 'v' && (buf[1] == 'n' || buf[1] == 't')))
		{
			fgets(buf, 256, pFile);
		}
		sscanf(buf, "v %lf %lf %lf", &x, &y, &z);
		pVertices[3 * i] = x;
		pVertices[3 * i + 1] = y;
		pVertices[3 * i + 2] = z;
	}

	/* read index data and average face normals */
	for (i = 0; i < iNumIndices; i += 3)
	{
		fgets(buf, 256, pFile);
		while (buf[0] != 'f')
		{
			fgets(buf, 256, pFile);
		}
		if (sscanf(buf, "f %d %d %d", pIndices + i, pIndices + i + 1, pIndices + i + 2) != 3)
		{
			int dmy[3];
			if (sscanf(buf, "f %d//%d %d//%d %d//%d", pIndices + i, dmy, pIndices + i + 1, dmy + 1, pIndices + i + 2, dmy + 2) != 6)
			{
				int dmy2[3];
				if (sscanf(buf, "f %d/%d/%d %d/%d/%d %d/%d/%d", pIndices + i, dmy2, dmy, pIndices + i + 1, dmy2 + 1, dmy + 1, pIndices + i + 2, dmy2 + 2, dmy + 2) != 9)
				{
					fprintf(stderr, "no support format!!\n");
					return;
				}
			}
		}
		pIndices[i] += -1;
		pIndices[i + 1] += -1;
		pIndices[i + 2] += -1;
		i0 = 3 * pIndices[i];
		i1 = 3 * pIndices[i + 1];
		i2 = 3 * pIndices[i + 2];
		VEC3_SUB(v1, pVertices + i1, pVertices + i0);
		VEC3_SUB(v2, pVertices + i2, pVertices + i0);
		VEC3_CROSS(n, v1, v2);
		VEC3_ADD(pNormals + i0, pNormals + i0, n);
		VEC3_ADD(pNormals + i1, pNormals + i1, n);
		VEC3_ADD(pNormals + i2, pNormals + i2, n);
	}
	fclose(pFile);

	/* normalize vertex normals */
	for (i = 0; i < iNumVertices; i += 3)
		VEC3_NORMALIZE(pNormals + i);

	printf("ok!\n");
	/* clean up */

	fp = fopen("tmp.off", "w");

	if (fp) fprintf(fp, "OFF\n");
	fprintf(fp, "%d %d 0\n", iNumVertices / 3, iNumIndices / 3);
	for (int i = 0; i < iNumVertices / 3; i++)
	{
		fprintf(fp, "%.10f %.10f %.10f\n", pVertices[3 * i], pVertices[3 * i + 1], pVertices[3 * i + 2]);
	}
	for (int i = 0; i < iNumIndices / 3; i++)
	{
		fprintf(fp, "3 %d %d %d\n", pIndices[3 * i], pIndices[3 * i + 1], pIndices[3 * i + 2]);
	}
	if (fp) fclose(fp);

	free(pVertices);
	free(pNormals);
	free(pIndices);

}
void load_off_to_obj(const char *filename, float rgb[3])
{
	unsigned int uiMesh;
	double *pVertices, *pNormals;
	unsigned int *pIndices;
	int i, j, f, iNumVertices, iNumIndices;
	unsigned int i0, i1, i2;
	double v1[3], v2[3], n[3];
	char buf[256];

	char tmpfile[256];
	char drive_[_MAX_DRIVE];	// ドライブ名
	char dir_[_MAX_DIR];		// ディレクトリ名
	char fname_[_MAX_FNAME];	// ファイル名
	char ext_[_MAX_EXT];		// 拡張子

	_splitpath(filename, drive_, dir_, fname_, ext_);

	FILE* fp = fopen(filename, "r");
	if (!fp)
	{
		fprintf(stderr, "cannot open file![%s]\n", filename);
		return;
	}


	fgets(buf, 256, fp);
	fgets(buf, 256, fp);
	int dmy;
	sscanf(buf, "%d %d %d\n", &iNumVertices, &iNumIndices, &dmy);

	pVertices = (double*)malloc(3*iNumVertices * sizeof(double));
	pNormals = (double*)calloc(3*iNumVertices, sizeof(double));
	pIndices = (unsigned int*)malloc(3*iNumIndices * sizeof(unsigned int));


	for (int i = 0; i < iNumVertices; i++)
	{
		fgets(buf, 256, fp);
		sscanf(buf, "%lf %lf %lf\n", &(pVertices[3 * i]), &(pVertices[3 * i + 1]), &(pVertices[3 * i + 2]));
	}
	for (int i = 0; i < iNumIndices; i++)
	{
		fgets(buf, 256, fp);
		sscanf(buf, "%d %d %d %d\n", &dmy, &(pIndices[3 * i]), &(pIndices[3 * i + 1]), &(pIndices[3 * i + 2]));
	}
	if (fp) fclose(fp);
	printf("vertex %d face %d\n", iNumVertices, iNumIndices);


	sprintf(tmpfile, "%s%s%s%s", drive_, dir_, fname_, ".obj");
	fp = fopen(tmpfile, "w");

	sprintf(tmpfile, "%s%s%s%s", drive_, dir_, fname_, ".mtl");
	fprintf(fp, "mtllib %s\n", tmpfile);
	for (i = 0; i < iNumVertices; ++i)
	{
		fprintf(fp, "v %.10f %.10f %.10f %f %f %f\n", 
			pVertices[3 * i], pVertices[3 * i + 1], pVertices[3 * i + 2],
			rgb[0], rgb[1], rgb[2]);
	}
	fprintf(fp, "usemtl mat\n");

	/* read index data and average face normals */
	for (i = 0; i < iNumIndices; ++i)
	{
		fprintf(fp, "f %d %d %d\n", pIndices[3 * i] + 1, pIndices[3 * i + 1] + 1, pIndices[3 * i + 2] + 1);
	}
	fclose(fp);

	sprintf(tmpfile, "%s%s%s%s", drive_, dir_, fname_, ".mtl");
	fp = fopen(tmpfile, "w");

	fprintf(fp, "newmtl mat\n");
	fprintf(fp, "Ka 0.25000 0.25000 0.25000\n");
	fprintf(fp, "Kd %f %f %f\n", rgb[0], rgb[1], rgb[2]);
	fprintf(fp, "Ks 0.25000 0.25000 0.25000\n");
	fprintf(fp, "Ns 5\n");
	fclose(fp);
	printf("ok!\n");

	/* clean up */
	free(pVertices);
	free(pNormals);
	free(pIndices);

}

int main(int argc, char** argv)
{
	char tmpfile[256];
	FILE* fp;
	char drive_[_MAX_DRIVE];	// ドライブ名
	char dir_[_MAX_DIR];		// ディレクトリ名
	char fname_[_MAX_FNAME];	// ファイル名
	char ext_[_MAX_EXT];		// 拡張子

	int 	number_of_clusters = 5;
	double smoothing_lambda = 0.26;

	for (int i = 2; i < argc-1; i++)
	{
		if (strcmp("-clusters", argv[i]) == 0)
		{
			number_of_clusters = atoi(argv[i + 1]);
			continue;
		}
		if (strcmp("-lambda", argv[i]) == 0)
		{
			smoothing_lambda = atof(argv[i + 1]);
			continue;
		}
	}
	SM mesh;
	if (argc >= 2) 
	{
		printf("input=[%s]\n", argv[1]);
		_splitpath(argv[1], drive_, dir_, fname_, ext_);

		if (std::string(ext_) == ".obj" || std::string(ext_) == ".OBJ")
		{
			load_obj_to_off(argv[1]);
			std::ifstream input("tmp.off");
			input >> mesh;
		}
		if (std::string(ext_) == ".off" || std::string(ext_) == ".OFF")
		{
			std::ifstream input(argv[1]);
			input >> mesh;
		}
	}
	else
	{
		fprintf(stderr, "auto_segmentation meshfilename [-clusters number_of_clusters] [-lambda smoothing_lambda]\n");
		return -1;
	}
	//else {
	//	std::ifstream cactus("data/cactus.off");
	//	cactus >> mesh;
	//}
	printf("number_of_faces:%d\n", mesh.number_of_faces());
	printf("number_of_vertices:%d\n", mesh.number_of_vertices());

	typedef SM::Property_map<face_descriptor, double> Facet_double_map;
	Facet_double_map sdf_property_map;

	sdf_property_map = mesh.add_property_map<face_descriptor, double>("f:sdf").first;

	CGAL::sdf_values(mesh, sdf_property_map);

	// create a property-map for segment-ids
	typedef SM::Property_map<face_descriptor, std::size_t> Facet_int_map;
	Facet_int_map segment_property_map = mesh.add_property_map<face_descriptor, std::size_t>("f:sid").first;;

	// segment the mesh using default parameters for number of levels, and smoothing lambda
	// Any other scalar values can be used instead of using SDF values computed using the CGAL function
	std::size_t number_of_segments = CGAL::segmentation_from_sdf_values(mesh, sdf_property_map, segment_property_map, number_of_clusters, smoothing_lambda);

	typedef CGAL::Face_filtered_graph<SM> Filtered_graph;
	//print area of each segment and then put it in a Mesh and print it in an OFF file
	Filtered_graph segment_mesh(mesh, 0, segment_property_map);

	printf("number_of_segments:%d\n", number_of_segments);

	std::mt19937 mt;
	std::uniform_real_distribution<> rand(0.2, 1.0);

	for (std::size_t id = 0; id < number_of_segments; ++id)
	{
		char tmpfile[256];

		if (id > 0)
			segment_mesh.set_selected_faces(id, segment_property_map);
		std::cout << "Segment " << id << "'s area is : " << CGAL::Polygon_mesh_processing::area(segment_mesh) << std::endl;
		SM out;
		CGAL::copy_face_graph(segment_mesh, out);
		std::ostringstream oss;
		sprintf(tmpfile, "%s%s%s_Segment_%d", drive_, dir_, fname_, id);
		//oss << tmpfile  << id << ".off";
		//oss << tmpfile << id;
		oss << tmpfile;
		std::ofstream os(oss.str().data());
		os << out;
		os.close();

		float rgb[3];

		rgb[0] = rand(mt);
		rgb[1] = rand(mt);
		rgb[2] = rand(mt);
		load_off_to_obj(tmpfile, rgb);
	}

}

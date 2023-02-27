#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <list>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

typedef Polyhedron::Facet_const_iterator Facet_iterator;
typedef Polyhedron::Vertex_const_iterator Vertex_iterator;
typedef Polyhedron::Halfedge_const_iterator Halfedge_iterator;
typedef Polyhedron::Halfedge_around_facet_const_circulator Halfedge_facet_circulator;

typedef std::map<Polyhedron::Facet_const_handle, double> Facet_double_map;
typedef std::map<Polyhedron::Facet_const_handle, int> Facet_int_map;

#define NBR_CLASSES 5

/// @brief map all the values from [min, max] to [0, 1]
/// @param facetMap non-const reference to the map (it is an in/out parameter)

struct Intervalle
{
	double borne_inf;
	double borne_sup;
};

// structure definie pour la question 3
struct Couleur
{
	double red;
	double green;
	double blue;
};

/*Retoune un vecteur de couleurs aléatoires*/
std::vector<Couleur> init_vecteur_couleurs(int nbr_classes)
{
	std::vector<Couleur> output;
	Couleur color = {0.0, 0.0, 0.0};
	int i = 0;
	do
	{
		color.red = ((double)rand() / (RAND_MAX));
		color.green = ((double)rand() / (RAND_MAX));
		color.blue = ((double)rand() / (RAND_MAX));
		output.push_back(color);
		// std::cout << "Color RGB is :" << color.red << " " << color.green << " " << color.blue << std::endl;
		i++;

	} while (i < nbr_classes + 1);
	std::cout << "fin de init_vecteur_couleurs" << std::endl;
	return output;
}

int nombre_facettes(Polyhedron &mesh)
{
	unsigned int nbFaces = 0;
	for (Facet_iterator i = mesh.facets_begin(); i != mesh.facets_end(); ++i)
	{
		++nbFaces;
	}
	return nbFaces;
}

void normalizeMap(Facet_double_map &facetMap)
{
	double maxValue = facetMap.begin()->second;
	double minValue = facetMap.begin()->second;

	// look for min and max value in the map
	for (const auto &elem : facetMap)
	{
		if (elem.second > maxValue)
		{
			maxValue = elem.second;
		}
		if (elem.second < minValue)
		{
			minValue = elem.second;
		}
	}

	for (auto &elem : facetMap)
	{
		elem.second -= minValue;
		elem.second /= (maxValue - minValue);
	}
}

/// @brief Generate in a .off file a colored mesh according to a value map (green to red shades)
/// @param mesh the input mesh
/// @param facetMap map of values between 0 and 1 (see "normalize()") for each facet of mesh
/// @param filePath path to the colored .off file to be generated
void writeOFFfromValueMap(const Polyhedron &mesh, const Facet_double_map &facetMap, std::string filePath)
{
	std::ofstream in_myfile;
	in_myfile.open(filePath);

	CGAL::set_ascii_mode(in_myfile);

	in_myfile << "COFF" << std::endl // "COFF" makes the file support color informations
			  << mesh.size_of_vertices() << ' '
			  << mesh.size_of_facets() << " 0" << std::endl;
	// nb of vertices, faces and edges (the latter is optional, thus 0)

	std::copy(mesh.points_begin(), mesh.points_end(),
			  std::ostream_iterator<Kernel::Point_3>(in_myfile, "\n"));

	for (Facet_iterator i = mesh.facets_begin(); i != mesh.facets_end(); ++i)
	{
		Halfedge_facet_circulator j = i->facet_begin();

		CGAL_assertion(CGAL::circulator_size(j) >= 3);

		in_myfile << CGAL::circulator_size(j) << ' ';
		do
		{
			in_myfile << ' ' << std::distance(mesh.vertices_begin(), j->vertex());

		} while (++j != i->facet_begin());

		in_myfile << std::setprecision(5) << std::fixed; // set the format of floats to X.XXXXX

		auto redValue = 1 - facetMap.at(i); // low values will be closer to red
		auto greenValue = facetMap.at(i);	// high values will be closer to green
		auto blueValue = 0.0;

		in_myfile << " " << redValue << " " << greenValue << " " << blueValue;

		in_myfile << std::endl;
	}

	in_myfile.close();

	std::cout << "Le résultat a été exporté dans " << filePath << " !" << std::endl;
}

Facet_double_map computePerimMap(const Polyhedron &mesh)
{
	Facet_double_map out;

	for (Facet_iterator i = mesh.facets_begin(); i != mesh.facets_end(); ++i)
	{
		double current_perimeter = 0.;
		Halfedge_facet_circulator j = i->facet_begin();
		do
		{
			current_perimeter += std::sqrt(CGAL::squared_distance(j->vertex()->point(), j->opposite()->vertex()->point()));
		} while (++j != i->facet_begin());

		std::cout << "perim(" << std::distance(mesh.facets_begin(), i) << ")=" << current_perimeter << std::endl;

		out[i] = current_perimeter;
	}

	return out;
}

Facet_double_map computeAreaMap(const Polyhedron &mesh)
{
	Facet_double_map out;

	for (Facet_iterator i = mesh.facets_begin(); i != mesh.facets_end(); ++i)
	{

		Halfedge_facet_circulator j = i->facet_begin();

		Polyhedron::Vertex_const_handle firstVertex = j->vertex();

		double current_area = 0;
		// a facet is not necessarily a triangle, so we decompose one facet into multiple triangles,
		// and sum up all their areas. Only works for convex faces.
		// (illustration: http://mathbitsnotebook.com/JuniorMath/Polygons/polygons3g.jpg)
		do
		{
			current_area += CGAL::squared_area(firstVertex->point(), j->vertex()->point(), j->opposite()->vertex()->point());
		} while (++j != i->facet_begin());

		// std::cout << "area(" << std::distance(mesh.facets_begin(), i) << ")=" << current_area << std::endl;

		out[i] = current_area;
	}

	return out;
}

// Question 2 du TP
Facet_int_map simpleThreshold(Polyhedron &mesh, Facet_double_map &values)
{

	Facet_int_map out;

	// calcul du seuil
	double seuil = 0.0;
	double somme = 0.0;
	int size = values.size();
	for (auto it = values.begin(); it != values.end(); ++it)
	{
		somme = somme + it->second;
		std::cout << "valeur : " << it->second << std::endl;
	}
	seuil = somme / size;
	std::cout << "seuil : " << seuil << std::endl;

	// attribution du seuillage
	for (auto it = values.begin(); it != values.end(); ++it)
	{
		if (it->second < seuil)
		{
			out[it->first] = 0;
		}
		else
		{
			out[it->first] = 1;
		}
	}
	return out;
}

// overload de la fonction pour un Facet_int_map
void writeOFFfromValueMap(const Polyhedron &mesh, const Facet_int_map &facetMap, std::string filePath, int nbr_classes)
{
	auto classe_max = std::max_element(facetMap.begin(), facetMap.end(),
									   [](const std::pair<Polyhedron::Facet_const_handle, int> &it1, const std::pair<Polyhedron::Facet_const_handle, int> &it2)
									   {
										   return it1.second < it2.second;
									   });

	std::vector<Couleur> color_vect = init_vecteur_couleurs(classe_max->second);
	//?auto color_vect_iterator = color_vect.begin();
	std::ofstream in_myfile;
	in_myfile.open(filePath);

	CGAL::set_ascii_mode(in_myfile);

	in_myfile << "COFF" << std::endl // "COFF" makes the file support color informations
			  << mesh.size_of_vertices() << ' '
			  << mesh.size_of_facets() << " 0" << std::endl;
	// nb of vertices, faces and edges (the latter is optional, thus 0)
	// std::cout << "deb" << std::endl;
	std::copy(mesh.points_begin(), mesh.points_end(),
			  std::ostream_iterator<Kernel::Point_3>(in_myfile, "\n"));

	for (Facet_iterator i = mesh.facets_begin(); i != mesh.facets_end(); ++i)
	{
		// std::cout << " started the ist for loop " << std::endl;
		Halfedge_facet_circulator j = i->facet_begin();

		CGAL_assertion(CGAL::circulator_size(j) >= 3);

		in_myfile << CGAL::circulator_size(j) << ' ';
		do
		{
			// std::cout << "do while" << std::endl;
			in_myfile << ' ' << std::distance(mesh.vertices_begin(), j->vertex());

		} while (++j != i->facet_begin());

		// std::cout << "apres do while" << std::endl;
		in_myfile << std::setprecision(6) << std::fixed; // set the format of floats to X.XXXXX

		// auto redValue = (double)1 - facetMap.at(i); // low values will be closer to red
		// std::cout << "red" << redValue << std::endl;
		// auto greenValue = (double)facetMap.at(i); // high values will be closer to green
		// std::cout << "green" << greenValue << std::endl;
		// auto blueValue = 0.0;

		/*****************************Partie ajoutée****************************************/

		auto index = facetMap.at(i); // valeurs entre 0 et nbr_classes+1
		// pour éviter un outofbound, le vecteur des couleurs contient pour l'instant nbr+1 couleurs
		// TODO : fix the problem!!

		std::cout << "facetMap.at(i) : " << facetMap.at(i) << std::endl;

		double redValue = color_vect.at(index).red;
		double greenValue = color_vect.at(index).green;
		double blueValue = color_vect.at(index).blue;

		// std::cout << "red : " << redValue << std::endl;
		// std::cout << "green : " << greenValue << std::endl;
		// std::cout << "blue : " << blueValue << std::endl;

		// std::cout << "Finished the ist  for loop " << std::endl;

		/*********************************************************************/

		in_myfile << " " << redValue << " " << greenValue << " " << blueValue;
		in_myfile << std::endl;
	}

	in_myfile.close();

	std::cout << "Le résultat a été exporté dans " << filePath << " !" << std::endl;
}

/**********************************************************/
// On retourne un tableau d'intervalles pour les classes auxquelles appartiennent les faces
std::vector<Intervalle> vecteur_intervalles(Facet_double_map &values, int nbr_classes)
{

	std::vector<Intervalle> classes;

	int i = 0;
	double pas = 0.0;
	double size = values.size();
	// double somme = std::accumulate(values.begin(), values.end(), 0,
	// 							   [](const double prev_sum, const std::pair<Polyhedron::Facet_const_handle, double> &entry)
	// 							   {
	// 								   return prev_sum + entry.second;
	// 							   });
	auto max_map = std::max_element(values.begin(), values.end(),
									[](const std::pair<Polyhedron::Facet_const_handle, double> &it1, const std::pair<Polyhedron::Facet_const_handle, double> &it2)
									{
										return it1.second < it2.second;
									});

	auto min_map = std::min_element(values.begin(), values.end(),
									[](const std::pair<Polyhedron::Facet_const_handle, double> &it1, const std::pair<Polyhedron::Facet_const_handle, double> &it2)
									{
										return it1.second < it2.second;
									});

	std::cout << "max_map : " << max_map->second << std::endl;
	std::cout << "min_map : " << min_map->second << std::endl;

	pas = (max_map->second - min_map->second) / nbr_classes;
	std::cout << "pas : " << pas << std::endl;

	// pas = max-min/nb_classes
	Intervalle interv = {min_map->second, min_map->second + pas};

	do
	{
		classes.push_back(interv);
		// std::cout << "intervalle actuel : [ " << interv.borne_inf
		// 		  << " , " << interv.borne_sup << "]" << std::endl;
		interv.borne_inf = interv.borne_sup;
		interv.borne_sup = interv.borne_inf + pas;
		i++;
	} while (i < nbr_classes);

	std::cout << "finished creating interval vector" << std::endl;

	return classes;
}

/*Retoune un itérateur sur la classe à laquelle coreespond la paire donnée*/
int chercher_classe(const std::pair<Polyhedron::Facet_const_handle, double> &paire, std::vector<Intervalle> &classes)
{
	int indice_classe = -1;
	auto politique = [paire](const Intervalle &interv)
	{ return (paire.second >= interv.borne_inf & paire.second < interv.borne_sup); };
	auto it = std::find_if(classes.begin(), classes.end(), politique);
	indice_classe = it - classes.begin();
	return indice_classe;
}

/// question3
Facet_int_map complexThreshold(Polyhedron &mesh, Facet_double_map &values, int nb_classes)
{
	std::ofstream myfile;
	myfile.open("out_segmentation_1.txt");
	myfile << "Writing this to a file.\n";

	Facet_int_map out;

	std::vector<Intervalle> classes = vecteur_intervalles(values, nb_classes);
	int indice_classe = -1;

	for (auto it = values.begin(); it != values.end(); ++it)
	{
		indice_classe = chercher_classe(*it, classes);
		out[it->first] = indice_classe;
		// std::cout << "out[it->first]" << out[it->first] << std::endl;
		//  out[it->second] = indice_classe; //error : no match for ‘operator[]’
		//  std::cout << "First : " << out[it->first] << std::endl;

		myfile << "out[it->first] : " << out[it->first] << std::endl;
	}

	return out;
}

// Facet_int_map segmentationParCC(Polyhedron &mesh, const Facet_int_map &segmentation)
// {
// 	Facet_int_map out;
// 	int nb_facettes = nombre_facettes(mesh);
// 	bool facette_visitee[nb_facettes];

// 	for (int i = 0; i < nb_facettes; i++)
// 	{
// 		facette_visitee[i] = 0;
// 	}

// 	//! Penser à le définir direct dans le main en var globale si besoin
// 	// std::vector<Intervalle> classes = vecteur_intervalles(values, nb_classes);
// 	int indice_i = 0;
// 	// parcourir les facettes
// 	for (Facet_iterator i = mesh.facets_begin(); i != mesh.facets_end(); ++i)
// 	{

// 		// auto indice_i = i - mesh.facets_begin();

// 		if (facette_visitee[indice_i] == 0) // facette non visitee
// 		{
// 			Halfedge_facet_circulator j = i->facet_begin();
// 			int indice_j = 0;
// 			do
// 			{

// 				// facette adjacente
// 				j->opposite()->facet();

// 				// si la facette a la même classe

// 				if (segmentation.at(i) == segmentation.at(j))
// 				{
// 					// appel récursif
// 				}
// 				else
// 				{
// 					// attribuer un nouvel id
// 				}

// 				indice_j++;

// 			} while (++j != i->facet_begin());
// 		}

// 		facette_visitee[indice_i] = 1; // on marque la facette comme visitée
// 		indice_i++;
// 	}

// 	return out;
// }

// retorne une structure des facettes adjascentes
std::vector<Polyhedron::Facet_const_handle> getFacettesAdjascentes(Facet_iterator i)
{

	std::vector<Polyhedron::Facet_const_handle> out;

	Halfedge_facet_circulator j = i->facet_begin();

	do
	{
		out.push_back(j->opposite()->facet());
	} while (++j != i->facet_begin());
	std::cout << "nb facette adj" << out.size() << std::endl;

	return out;
}

void parcours_profond(Polyhedron &mesh, Facet_int_map &segmentation, Facet_iterator i, bool *tableau_parcours, int classe, Facet_int_map &out)
{
	Facet_iterator face_deb = mesh.facets_begin();

	std::vector<Polyhedron::Facet_const_handle> facettes_adjascentes = getFacettesAdjascentes(i);

	std::cout << "std::distance(face_deb, i)" << std::distance(face_deb, i) << std::endl;

	if (tableau_parcours[std::distance(face_deb, i)]) // visitee
	{
		std::cout << "Facette deja visitee dans ParcoursProfond" << std::endl;

		return;
	}
	else
	{
		out[i] = classe;
		std::cout << " classe : " << classe << std::endl;
		tableau_parcours[std::distance(face_deb, i)] = true;
		for (auto val : facettes_adjascentes)
		{
			if ((segmentation[i] == segmentation[val]) && (!tableau_parcours[std::distance(face_deb, val)])) // sinon boucle infinie
			{
				out[val] = classe;
				parcours_profond(mesh, segmentation, val, tableau_parcours, classe, out);
				// tableau_parcours[std::distance(face_deb, val)] = true;
			}
		}
	}
}

Facet_int_map segmentationParCC(Polyhedron &mesh, Facet_int_map &segmentation)
{

	std::ofstream myfile;
	myfile.open("out_segmentation_2.txt");
	myfile << "Writing this to a file.\n";

	Facet_int_map out;

	Facet_iterator face_deb = mesh.facets_begin();

	bool tableau_parcours[mesh.size_of_facets()] = {false};
	int classe = 0;

	for (Facet_iterator i = mesh.facets_begin(); i != mesh.facets_end(); ++i)
	{

		// std::cout << "for std::distance(face_deb, i) : " << std::distance(face_deb, i)
		// 		  << " tableau_parcours[std::distance(face_deb, i)] : " << tableau_parcours[std::distance(face_deb, i)]
		// 		  << std::endl;

		if (tableau_parcours[std::distance(face_deb, i)]) //  facette visitee
		{
			std::cout << "Facette deja visitee dans Segmentation cc" << std::endl;
			// break;
		}
		else // facette non visitee
		{
			// tableau_parcours[std::distance(face_deb, i)] = true;
			parcours_profond(mesh, segmentation, i, tableau_parcours, classe, out);
			// out[i] = classe;
			++classe;
		}

		// current_class = segmentation[i];
		// std::cout << "out[i] : " << out[i] << std::endl;

		myfile << "out[i] : " << out[i] << std::endl;
	}
	std::cout << "ajout_classe : " << classe << std::endl;

	return out;
}

int main(int argc, char *argv[])
{
	if (argc < 2)
	{
		std::cerr << "Il manque un paramètre au programme. Veuillez lui donner en entrée un nom de fichier au format off." << std::endl;
		return 1;
	}

	Polyhedron mesh;

	std::ifstream input(argv[1]);

	if (!input || !(input >> mesh) || mesh.is_empty())
	{
		std::cerr << "Le fichier donné n'est pas un fichier off valide." << std::endl;
		return 1;
	}

	// auto mapPerim = computePerimMap(mesh);

	auto mapArea = computeAreaMap(mesh);

	normalizeMap(mapArea);

	// writeOFFfromValueMap(mesh, mapPerim, argc >= 3 ? argv[2] : "result.off");

	// auto map_test = simpleThreshold(mesh, mapPerim);

	// writeOFFfromValueMap(mesh, map_test, argc >= 3 ? argv[2] : "question2.off");

	// auto map_test_2 = simpleThreshold(mesh, mapPerim, 3);

	// auto testing = vecteur_intervalles(mapPerim, 10.0);

	auto map_test_3 = complexThreshold(mesh, mapArea, NBR_CLASSES);
	writeOFFfromValueMap(mesh, map_test_3, argc >= 3 ? argv[2] : "question3.off", NBR_CLASSES);

	// // std::cout << "Size de la map 3 avant segmentationParCC : " << map_test_3.size() << std::endl;
	auto map_test_4 = segmentationParCC(mesh, map_test_3);
	// // std::cout << "Size de la map 4 après segmentationParCC: " << map_test_4.size() << std::endl;

	writeOFFfromValueMap(mesh, map_test_4, argc >= 3 ? argv[2] : "question4.off", NBR_CLASSES);
	//    auto testing_color_vector = init_vecteur_couleurs(10);

	return 0;
}

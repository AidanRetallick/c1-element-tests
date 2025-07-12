//==============================================================================
//==============================================================================
// Test file to output the basis functions for each CurveableBellElements
//==============================================================================
//==============================================================================

#include <iostream>
#include "c1_foeppl_von_karman/foeppl_von_karman_curvable_bell_elements.h"

using namespace oomph;


// COPIED FROM MATTHIAS' DISC EXPLORER CODE
 //=========================================================================
 /// Steady, straight 1D line in 2D space connecting two specified points.
 /// First point reached for zeta = 0; last one for zeta = 1 
 //=========================================================================
class TwoDStraightLineFromTwoPoints : public GeomObject
{
public:
 
 /// Constructor: Pass left and right point.
 TwoDStraightLineFromTwoPoints(const Vector<double>& left,
                               const Vector<double>& right) 
  : GeomObject(1, 2)
  {
#ifdef PARANOID
   if (left.size() != 2)
    {
     std::ostringstream error_message;
     error_message << "left point should have size 2, not "
                   << left.size() << std::endl;
     throw OomphLibError(error_message.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
   if (right.size() != 2)
    {
     std::ostringstream error_message;
     error_message << "right point should have size 2, not "
                   << right.size() << std::endl;
     throw OomphLibError(error_message.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
#endif

   Left.resize(2);
   Left[0]=left[0];
   Left[1]=left[1];
   Right.resize(2);
   Right[0]=right[0];
   Right[1]=right[1];
   
  }
 
 /// Broken copy constructor
 TwoDStraightLineFromTwoPoints(const TwoDStraightLineFromTwoPoints& dummy) = delete;
 
 /// Broken assignment operator
 void operator=(const TwoDStraightLineFromTwoPoints&) = delete;
 
 /// Destructor
 ~TwoDStraightLineFromTwoPoints(){}
 
 /// Position Vector at Lagrangian coordinate zeta
 void position(const Vector<double>& zeta, Vector<double>& r) const
  {
   // Position Vector
   r[0] = Left[0]+zeta[0]*(Right[0]-Left[0]);
   r[1] = Left[1]+zeta[0]*(Right[1]-Left[1]);
  }
 
 
 /// Parametrised position on object: r(zeta). Evaluated at
 /// previous timestep. t=0: current time; t>0: previous
 /// timestep.
 void position(const unsigned& t,
               const Vector<double>& zeta,
               Vector<double>& r) const
  {
   // Position Vector
   r[0] = Left[0]+zeta[0]*(Right[0]-Left[0]);
   r[1] = Left[1]+zeta[0]*(Right[1]-Left[1]);
  }
 
 
 /// Derivative of position Vector w.r.t. to coordinates:
 /// \f$ \frac{dR_i}{d \zeta_\alpha}\f$ = drdzeta(alpha,i).
 /// Evaluated at current time.
 virtual void dposition(const Vector<double>& zeta,
                        DenseMatrix<double>& drdzeta) const
  {
   // Tangent vector
   drdzeta(0, 0) = Right[0]-Left[0];
   drdzeta(0, 1) = Right[1]-Left[1];
  }
 
 
 /// 2nd derivative of position Vector w.r.t. to coordinates:
 /// \f$ \frac{d^2R_i}{d \zeta_\alpha d \zeta_\beta}\f$ =
 /// ddrdzeta(alpha,beta,i). Evaluated at current time.
 virtual void d2position(const Vector<double>& zeta,
                         RankThreeTensor<double>& ddrdzeta) const
  {
   // Derivative of tangent vector
   ddrdzeta(0, 0, 0) = 0.0;
   ddrdzeta(0, 0, 1) = 0.0;
  }
 
 
 /// Posn Vector and its  1st & 2nd derivatives
 /// w.r.t. to coordinates:
 /// \f$ \frac{dR_i}{d \zeta_\alpha}\f$ = drdzeta(alpha,i).
 /// \f$ \frac{d^2R_i}{d \zeta_\alpha d \zeta_\beta}\f$ =
 /// ddrdzeta(alpha,beta,i).
 /// Evaluated at current time.
 virtual void d2position(const Vector<double>& zeta,
                         Vector<double>& r,
                         DenseMatrix<double>& drdzeta,
                         RankThreeTensor<double>& ddrdzeta) const
  {
   // Position Vector
   r[0] = Left[0]+zeta[0]*(Right[0]-Left[0]);
   r[1] = Left[1]+zeta[0]*(Right[1]-Left[1]);
   
   // Tangent vector
   drdzeta(0, 0) = Right[0]-Left[0];
   drdzeta(0, 1) = Right[1]-Left[1];
   
   // Derivative of tangent vector
   ddrdzeta(0, 0, 0) = 0.0;
   ddrdzeta(0, 0, 1) = 0.0;
  }
 
 
 /// How many items of Data does the shape of the object depend on?
 unsigned ngeom_data() const
  {
   return 0;
  }
 
 /// Return pointer to the j-th Data item that the object's
 /// shape depends on
 Data* geom_data_pt(const unsigned& j)
  {
   return 0;
  }

private:
 
 /// Left point 
 Vector<double> Left;

 /// Right point
 Vector<double> Right; 
 
};

int main()
{
  // Define vertices
  Vector<double> v0{0.0, 2.0};
  Vector<double> v1{0.0, 0.0};
  Vector<double> v2{1.0, 1.0};

  // Create the elements (suffixed based on boundary interpolation)
  FoepplVonKarmanC1CurvableBellElement<4> el_straight;
  FoepplVonKarmanC1CurvableBellElement<4> el_cubic;
  FoepplVonKarmanC1CurvableBellElement<4> el_quintic;
  // Array of element pointers
  // (named for what it is Matthias, not what it contains)
  FoepplVonKarmanC1CurvableBellElement<4>* el_pt_array[3] =
    {&el_straight, &el_cubic, &el_quintic};
  for(auto el_pt : el_pt_array)
  {
    Node* n0 = el_pt->construct_node(0);
    Node* n1 = el_pt->construct_node(1);
    Node* n2 = el_pt->construct_node(2);
    n0->x(0) = v0[0];
    n0->x(1) = v0[1];
    n1->x(0) = v1[0];
    n1->x(1) = v1[1];
    n2->x(0) = v2[0];
    n2->x(1) = v2[1];
  }

  // Create the straight curvilinear boundary object and upgrade the two
  // higher order elements
  TwoDStraightLineFromTwoPoints* straight_edge_pt =
    new TwoDStraightLineFromTwoPoints(v0, v1);

  TriangleMeshCurviLine* temp_straight_edge_pt =
    new TriangleMeshCurviLine(straight_edge_pt, 0.0, 1.0, 1, 0);

  C1CurviLine* curviline_straight_edge_pt =
    new C1CurviLine(temp_straight_edge_pt);

  C1PlateHelper::CurvedEdgeEnumeration edge_number =
    C1PlateHelper::CurvedEdgeEnumeration::two;
  el_cubic.upgrade_element_to_curved(edge_number, 0, 1, curviline_straight_edge_pt, 3);
  el_quintic.upgrade_element_to_curved(edge_number, 0, 1, curviline_straight_edge_pt, 5);


  //
  char reslt_dir[7] = "RESLT";
  // Number of each type
  unsigned n_node = 3;
  unsigned n_nodal_type = 6;
  // Loop over each element
  for(unsigned i_el = 0; i_el < 3; i_el++)
  {
    FoepplVonKarmanC1CurvableBellElement<4>* el_pt = el_pt_array[i_el];
    unsigned n_internal_type = el_pt->ninternal_basis_type_for_field(2);
    // Loop over each of the functions and document
    for (unsigned i = 0; i < n_node; i++)
    {
      for (unsigned j = 0; j < n_nodal_type; j++)
      {
          unsigned w_index = 2;
          char filename[100];
          std::sprintf(filename, "%s/element_%i_%i.dat", reslt_dir, n_internal_type, 6*i + j);
          std::ofstream outfile(filename);
          el_pt->node_pt(i)->set_value(w_index + j, 1.0);
          el_pt->full_output(outfile, 100);
          el_pt->node_pt(i)->set_value(w_index + j, 0.0);
          outfile.close();
      }
    }
    for (unsigned j = 0; j < n_internal_type; j++)
    {
        char filename[100];
        std::sprintf(filename, "%s/element_%i_%i.dat", reslt_dir, n_internal_type, 18 + j);
        std::ofstream outfile(filename);
        el_pt->internal_data_pt(2)->set_value(j, 1.0);
        el_pt->full_output(outfile, 100);
        el_pt->internal_data_pt(2)->set_value(j, 0.0);
        outfile.close();
    }
  }

  return 0;
}

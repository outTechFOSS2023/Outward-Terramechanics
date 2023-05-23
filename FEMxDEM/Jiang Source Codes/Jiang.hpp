// 2017 Nima Goudarzi ©  <nima.goudarzi1@northwestern.edu> 
// 
/*
=== A NOVEL THREE-DIMENSIONAL CONTACT MODEL FOR GRANULATES INCORPORATING ROLLING AND TWISTING RESISTANCES ===

This is implementation of a 3D contact model incorporating rolling and twisting resistances at inter-particle contact, 
to simulate the mechanical behavior of particulates.
The DMT formulation is also considered (for adhesive particles, rigid and small bodies).

*/


#pragma once

#include <pkg/dem/FrictPhys.hpp>
#include <pkg/common/ElastMat.hpp>
#include <pkg/common/Dispatching.hpp>
#include <pkg/dem/ScGeom.hpp>
#include <pkg/common/PeriodicEngines.hpp>
#include <pkg/common/NormShearPhys.hpp>
#include <pkg/common/MatchMaker.hpp>

#include <boost/tuple/tuple.hpp>
#include <lib/base/openmp-accu.hpp>


/************************************************************/
/************************ JiangMat **************************/
/************************************************************/
class JiangMat : public FrictMat
{
	public :
		virtual ~JiangMat () {};
/// Serialization
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(JiangMat,FrictMat,"",
		((Real,beta,0.0,,"Dimensionless shape parameter linking the contact radius $\\R_bar$ to the common radius r"))
		((Real,xIc,2.1,,"Local crushing parameter describing the effects of local asperity crushing "))
		,
		createIndex();
		);
/// Indexable
	REGISTER_CLASS_INDEX(JiangMat,FrictMat);
};

REGISTER_SERIALIZABLE(JiangMat);



/************************************************************/
/************************* JiangPhys ************************/
/************************************************************/
class JiangPhys: public FrictPhys{
	public:
	virtual ~JiangPhys() {};
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(JiangPhys,FrictPhys,"Representation of an interaction of the Jiang type.",
			((Real,kr,0.0,,"Rotational stiffness"))
			((Real,ktw,0.0,,"Rotational stiffness"))
			((Real,XIC,2.1,,"Local crushing parameter describing the effects of local asperity crushing"))
			((Real,R_bar,0.0,,"Contact radius as the product of common radius,r and beta"))
			
			((Vector3r,normalViscous,Vector3r::Zero(),,"Normal viscous component"))
			((Vector3r,shearViscous,Vector3r::Zero(),,"Shear viscous component"))
			((Vector3r,rollViscous,Vector3r::Zero(),,"Rolling viscous component"))
			((Vector3r,twistViscous,Vector3r::Zero(),,"Twisting viscous component"))
			((Vector3r,shearElastic,Vector3r::Zero(),,"Total elastic shear force"))
			((Vector3r,usElastic,Vector3r::Zero(),,"Total elastic shear displacement (only elastic part)"))
			((Vector3r,usTotal,Vector3r::Zero(),,"Total elastic shear displacement (elastic+plastic part)"))
			((Vector3r,rollMomentElastic,Vector3r::Zero(),,"Total elastic rolling moment"))
			((Vector3r,tetaRollElastic,Vector3r::Zero(),,"Total elastic rolling rotation (only elastic part)"))
			((Vector3r,tetaRollTotal,Vector3r::Zero(),,"Total elastic rolling rotation (elastic+plastic part)"))
			((Vector3r,twistMomentElastic,Vector3r::Zero(),,"Total elastic twisting moment"))
			((Vector3r,tetaTwistElastic,Vector3r::Zero(),,"Total elastic twisting rotation (only elastic part)"))
			((Vector3r,tetaTwistTotal,Vector3r::Zero(),,"Total elastic twisting rotation (elastic+plastic part)"))		
			//((Vector3r,dThetaR,Vector3r::Zero(),,"Incremental rolling vector"))
			((Vector3r,rollMoment,Vector3r::Zero(),,"Artificial rolling moment to provide rolling resistance in order to account for some degree of interlocking between particles"))
			((Vector3r,twistMoment,Vector3r::Zero(),,"Artificial twisting moment to provide rolling resistance in order to account for some degree of interlocking between particles"))
			//((Vector3r,prevNormal,Vector3r::Zero(),,"Save previous contact normal to compute relative rotation"))

			((Real,radius,NaN,,"Contact radius (only computed with :yref:`Law2_ScGeom_JiangPhys_Jiang::calcEnergy`)"))

			//((Real,gamma,0.0,"Surface energy parameter [J/m^2] per each unit contact surface, to derive DMT formulation from HM"))
			((Real,adhesionForce,0.0,,"Force of adhesion as predicted by DMT"))
			((bool,isAdhesive,false,,"bool to identify if the contact is adhesive, that is to say if the contact force is attractive"))
			((bool,isSliding,false,,"check if the contact is sliding (useful to calculate the ratio of sliding contacts)"))
			((bool,therRollRotReached,false,,"check if the threshold rolling rotation has been reached)"))
			((bool,therTwistRotReached,false,,"check if the threshold twisting rotation has been reached)"))

			// Direct contact damping ratio when betan and/or betas is/are given
			((Real,betanDir,0.0,,"Normal Damping Ratio. Fraction of the viscous damping coefficient (normal direction) equal to $\\frac{c_{n}}{C_{n,crit}}$."))
			((Real,betasDir,0.0,,"Shear Damping Ratio. Fraction of the viscous damping coefficient (shear direction) equal to $\\frac{c_{s}}{C_{s,crit}}$."))

			// Indirect contact damping ratio when en and/or es is/are given
			((Real,betanIndir,0.0,,"Normal Damping Ratio calculated from $e_n$ as $\\beta_n=-(\\log e_n)/\\sqrt{\\pi^2+(\\log e_n)^2}$."))
			

			,
			createIndex());
	REGISTER_CLASS_INDEX(JiangPhys,FrictPhys);
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(JiangPhys);



/************************************************************/
/**************Ip2_JiangMat_JiangMat_JiangPhys **************/
/************************************************************/
class Ip2_JiangMat_JiangMat_JiangPhys: public IPhysFunctor{
	public :
	virtual void go(const shared_ptr<Material>& b1,	const shared_ptr<Material>& b2,	const shared_ptr<Interaction>& interaction);
	FUNCTOR2D(JiangMat,JiangMat);
	YADE_CLASS_BASE_DOC_ATTRS(
			Ip2_JiangMat_JiangMat_JiangPhys,IPhysFunctor,"Calculate some physical parameters needed to obtain \
the normal and shear stiffnesses.\n\n\
Viscous parameters can be specified either using coefficients of restitution ($e_n$, $e_s$) or viscous \
damping ratio ($\\beta_n$, $\\beta_s$). The following rules apply:\n#. If the $\\beta_n$ ($\\beta_s$) \
ratio is given, it is assigned to :yref:`JiangPhys.betanDir` directly (it is assumed that $e_n$ and $e_n$ are equal).\n#. \
If $e_n$ is given, :yref:`JiangPhys.betanIndir` is computed using $\\beta_n=-(\\log e_n)/\\sqrt{\\pi^2+(\\log e_n)^2}$. \
The same applies to $e_s$, :yref:`JiangPhys.betas`.\n#. It is an error (exception) to specify both $e_n$ \
and $\\beta_n$ ($e_s$ and $\\beta_s$).\n#. If neither $e_n$ nor $\\beta_n$ is given, zero value \
for :yref:`JiangPhys.betan` is used; there will be no viscous effects.\n#.If neither $e_s$ nor $\\beta_s$ \
is given, the value of :yref:`JiangPhys.betan` is used for :yref:`JiangPhys.betas` as well.\n\nThe \
$e_n$, $\\beta_n$, $e_s$, $\\beta_s$ are :yref:`MatchMaker` objects; they can be constructed from float \
values to always return constant value.\n\nSee :ysrc:`scripts/test/shots.py` for an example of specifying \
$e_n$ based on combination of parameters.",
			((Real,gamma,0.0,,"Surface energy parameter [J/m^2] per each unit contact surface, to derive DMT formulation from HM"))
			//((Real,krot,0.0,,"Rotational stiffness for moment contact law (since kr in Iphys is directly calculated from kn, we don't need krot"))
			//((Real,ktwist,0.0,,"Torsional stiffness for moment contact law (since ktw in Iphys is directly calculated from ks, we don't need ktwist"))
			((shared_ptr<MatchMaker>,en,,,"Normal coefficient of restitution $e_n$."))
			((shared_ptr<MatchMaker>,es,,,"Shear coefficient of restitution $e_s$."))
			((shared_ptr<MatchMaker>,betan,,,"Normal viscous damping ratio $\\beta_n$."))
			((shared_ptr<MatchMaker>,betas,,,"Shear viscous damping ratio $\\beta_s$."))
			((shared_ptr<MatchMaker>,frictAngle,,,"Instance of :yref:`MatchMaker` determining how to compute the friction angle of an interaction. If ``None``, minimum value is used."))
	);
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(Ip2_JiangMat_JiangMat_JiangPhys);


/**********************************************************/
/************** Law2_ScGeom_JiangPhys_Jiang ***************/
/**********************************************************/
class Law2_ScGeom_JiangPhys_Jiang: public LawFunctor{
	public:

		virtual bool go(shared_ptr<IGeom>& _geom, shared_ptr<IPhys>& _phys, Interaction* I);
		Real adhesionEnergy(); 	
		Real normElastEnergy();
		//there is no plastic dissipation associated with normal direction
		Real getnormDampDissip();
		Real getshearEnergy();
		Real getfrictionDissipation();
		Real getshearDampDissip();
		Real getrollEnergy();
		Real getrollingPlasticDissipation();
		Real getrollingDampDissip();
		Real gettwistEnergy();
		Real gettwistPlasticDissipation();
		Real gettwistDampDissip();
		Real contactsAdhesive();
		Real ratioSlidingContacts();

		FUNCTOR2D(ScGeom,JiangPhys);
		YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(Law2_ScGeom_JiangPhys_Jiang,LawFunctor,"Constitutive law for the Jiang formulation. It includes linear elasticity in the normal direction for two non-conforming elastic contact bodies. In the shear direction, as well, a linear relationship between shear force and tangential displacement up to a maxFs is provided. After this, the Mohr-Coulomb criterion is employed to establish the maximum friction force which can be developed at the contact. This model, moreover, employees the same principle in rolling and twisting direction where a linear relationship between rolling (or twisting) moment and rolling (or twisting) rotation up to a maxMr (or maxMt) is provided. After this, a Mohr-Coulomb like criterion is employed to establish the maximum rolling (or twisting) moment which can be transferred at the contact. Moreover, it is also possible to include both direct (through the definition of the parameters $\\beta_{n}$ and $\\beta_{s}$) and indirect (through the definition of the parameters ($e_n$ and $e_s$) viscous damping in all for directions as explained in the application (.cpp) file.",
			((bool,preventGranularRatcheting,true,,"bool to avoid granular ratcheting"))
			((bool,includeAdhesion,false,,"bool to include the adhesion force following the DMT formulation. If true, also the normal elastic energy takes into account the adhesion effect."))
			((bool,calcEnergy,false,,"bool to calculate energy terms (shear,rooling and twisting potential energy, energy dissipation due to normal, shear, rolling and twisting viscous damping and shear, plastic dissipation in shear, rolling and twisting directions)"))
			((bool,includeRollResistMoment,false,,"bool to consider rolling resistance moment.)"))
			((bool,includeTwistResistMoment,false,,"bool to consider twisting resistance moment.)"))			
			((bool,neverErase,false,,"Keep interactions even if particles go away from each other (only in case another constitutive law is in the scene, e.g. :yref:`Law2_ScGeom_CapillaryPhys_Capillarity`)"))

			((OpenMPAccumulator<Real>,normDampDissip,,Attr::noSave,"Energy dissipation due to normal viscous damping"))
			((OpenMPAccumulator<Real>,shearEnergy,,Attr::noSave,"Shear elastic potential energy"))
			((OpenMPAccumulator<Real>,frictionDissipation,,Attr::noSave,"Energy dissipation due to plastic state (sliding) in tangential direction"))		
			((OpenMPAccumulator<Real>,shearDampDissip,,Attr::noSave,"Energy dissipation due to shear viscous damping"))
			((OpenMPAccumulator<Real>,rollEnergy,,Attr::noSave,"Rolling elastic potential energy"))
			((OpenMPAccumulator<Real>,rollingPlasticDissipation,,Attr::noSave,"Energy dissipation due to plastic state in rolling direction"))		
			((OpenMPAccumulator<Real>,rollingDampDissip,,Attr::noSave,"Energy dissipation due to rolling viscous damping"))
			((OpenMPAccumulator<Real>,twistEnergy,,Attr::noSave,"Twisting elastic potential energy"))
			((OpenMPAccumulator<Real>,twistPlasticDissipation,,Attr::noSave,"Energy dissipation due to plastic state in twisting direction"))		
			((OpenMPAccumulator<Real>,twistDampDissip,,Attr::noSave,"Energy dissipation due to twisting viscous damping"))			
			, /* init */
			, /* ctor */
			, /* py */
			.def("contactsAdhesive",&Law2_ScGeom_JiangPhys_Jiang::contactsAdhesive,"Compute total number of adhesive contacts.")
			.def("ratioSlidingContacts",&Law2_ScGeom_JiangPhys_Jiang::ratioSlidingContacts,"Return the ratio between the number of contacts sliding to the total number at a given time.")
			.def("normElastEnergy",&Law2_ScGeom_JiangPhys_Jiang::normElastEnergy,"Compute normal elastic potential energy. It handles the DMT formulation if :yref:`Law2_ScGeom_JiangPhys_Jiang::includeAdhesion` is set to true.")	
	);
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(Law2_ScGeom_JiangPhys_Jiang);


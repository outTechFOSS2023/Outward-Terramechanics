// 2017 Nima Goudarzi ©  <nima.goudarzi1@northwestern.edu> 

#include"Jiang.hpp"
#include<pkg/dem/ScGeom.hpp>
#include<core/Omega.hpp>
#include<core/Scene.hpp>
#include<math.h>

YADE_PLUGIN(
	(JiangMat)
	(JiangPhys)
	(Ip2_JiangMat_JiangMat_JiangPhys)
	(Law2_ScGeom_JiangPhys_Jiang)

);

Real Law2_ScGeom_JiangPhys_Jiang::getnormDampDissip() {return (Real) normDampDissip;}
Real Law2_ScGeom_JiangPhys_Jiang::getshearEnergy() {return (Real) shearEnergy;}
Real Law2_ScGeom_JiangPhys_Jiang::getfrictionDissipation() {return (Real) frictionDissipation;}
Real Law2_ScGeom_JiangPhys_Jiang::getshearDampDissip() {return (Real) shearDampDissip;}
Real Law2_ScGeom_JiangPhys_Jiang::getrollEnergy() {return (Real) rollEnergy;}
Real Law2_ScGeom_JiangPhys_Jiang::getrollingPlasticDissipation() {return (Real) rollingPlasticDissipation;}
Real Law2_ScGeom_JiangPhys_Jiang::getrollingDampDissip() {return (Real) rollingDampDissip;}
Real Law2_ScGeom_JiangPhys_Jiang::gettwistEnergy() {return (Real) twistEnergy;}
Real Law2_ScGeom_JiangPhys_Jiang::gettwistPlasticDissipation() {return (Real) twistPlasticDissipation;}
Real Law2_ScGeom_JiangPhys_Jiang::gettwistDampDissip() {return (Real) twistDampDissip;}

/************************************************************/
/**************Ip2_JiangMat_JiangMat_JiangPhys **************/
/************************************************************/

CREATE_LOGGER(Ip2_JiangMat_JiangMat_JiangPhys);

void Ip2_JiangMat_JiangMat_JiangPhys::go(const shared_ptr<Material>& b1,const shared_ptr<Material>& b2, const shared_ptr<Interaction>& interaction){
	if(interaction->phys) return; // no updates of an already existing contact necessary
	shared_ptr<JiangPhys> contactPhysics(new JiangPhys());
	interaction->phys = contactPhysics;
	JiangMat* mat1 = YADE_CAST<JiangMat*>(b1.get());
	JiangMat* mat2 = YADE_CAST<JiangMat*>(b2.get());
	
	/* from interaction physics */
	Real Ea = mat1->young;
	Real Eb = mat2->young; //only Ea will be used in calculation of Knorm
	Real Va = mat1->poisson; //this is xi for calculation of Kshear
	Real Vb = mat2->poisson; //(may be needed later)
	Real fa = mat1->frictionAngle;
	Real fb = mat2->frictionAngle;
	Real Beta = mat1->beta;
	Real Xic = mat1->xIc; // Local crushing parameter, (see JiangMat and top of page 151 of Jiang's paper definition of xIc)


	/* from interaction geometry */
	ScGeom* scg = YADE_CAST<ScGeom*>(interaction->geom.get());		
	Real Da = scg->radius1; 
	Real Db = scg->radius2; 
	//Vector3r normal=scg->normal;        //The variable set but not used


	/* calculate stiffness coefficients */
	//Real Ga = Ea/(2*(1+Va));
	//Real Gb = Eb/(2*(1+Vb));
	//Real G = (Ga+Gb)/2; // average of shear modulus
	Real V = (Va+Vb)/2; // average of poisson's ratio
	Real E = Ea*Eb/((1.-pow(Va,2))*Eb+(1.-pow(Vb,2))*Ea); // Young modulus
	Real R = Da*Db/(Da+Db); // equivalent radius (will be used in calculation of DMT adhesion. This is from original HM but will be kept for now)
	Real r = (2*Da*Db)/(Da+Db); // common radius (equation 35 of Jiang's paper)
	Real Knorm = 2*r*E; // calculte normal stiffness (equation 36a of Jiang's paper)
	Real Kshear = Knorm/V; // calculate shear stiffness (equation 36b of Jiang's paper)
	Real r_Bar = r*Beta; // contact radius (see JiangMat and equation 34 of Jiang's paper for definition of beta)
	Real Kroll = 0.25*Knorm*pow(r_Bar,2); // calculate rolling stiffness (equation 10 of Jiang's paper).In this model, Kroll is not a free parameter (it depends on Knorm)
	Real ktwist = 0.5*Kshear*pow(r_Bar,2); // calculate twisting stiffness (equation 22 of Jiang's paper). In this model, ktwist is not a free parameter (it depends on Kshear)
	Real frictionAngle = (!frictAngle) ? std::min(fa,fb) : (*frictAngle)(mat1->id,mat2->id,mat1->frictionAngle,mat2->frictionAngle);// as in original HM
    //Real frictionAngle = std::min(fa,fb) // for simlicity, this can also be used

	Real Adhesion = 4.*M_PI*R*gamma; // calculate adhesion force as predicted by DMT theory (for now use original HM using R not r)

	/* pass values calculated from above to JiangPhys */
	contactPhysics->tangensOfFrictionAngle = std::tan(frictionAngle); 
	//contactPhysics->prevNormal = scg->normal; // used to compute relative rotation (this is required in calculation of dThetaR(see Iphys for definition of dThetaR))
	contactPhysics->kn = Knorm; // pass normal stiffness to Iphys 
	contactPhysics->ks = Kshear; // pass shear stiffness to Iphys 
	contactPhysics->kr = Kroll; // pass rolling stiffness to Iphys 
	contactPhysics->ktw = ktwist; // pass twisting stiffness to Iphys 
	contactPhysics->adhesionForce = Adhesion;
	
	contactPhysics->R_bar =r_Bar; // pss contact radius to Iphys
	contactPhysics->XIC =Xic; // pass local crushing parameter, (see JiangMat and page 151 of Jiang's paper for definition of xIc )

	/* compute viscous coefficients */
	if(en && betan) throw std::invalid_argument("Ip2_JiangMat_JiangMat_JiangPhys: only one of en, betan can be specified.");
	if(es && betas) throw std::invalid_argument("Ip2_JiangMat_JiangMat_JiangPhys: only one of es, betas can be specified.");

	// en or es specified, just compute betanIndir, otherwise betanIndir remains 0
	if(en || es){
		Real logE = log((*en)(mat1->id,mat2->id));
		contactPhysics->betanIndir = -logE/sqrt(pow(logE,2)+pow(M_PI,2));
	}
	
	// betan specified, use that value directly; otherwise give zero
	else{	
		contactPhysics->betanDir=betan ? (*betan)(mat1->id,mat2->id) : 0; 
		contactPhysics->betasDir=betas ? (*betas)(mat1->id,mat2->id) : contactPhysics->betanDir;
	}
}

/*****************************************************************************************/
/* FUNCTION TO COUNT THE NUMBER OF ADHESIVE CONTACTS IN THE SIMULATION AT EACH TIME STEP */
/*****************************************************************************************/

Real Law2_ScGeom_JiangPhys_Jiang::contactsAdhesive() // It is returning something rather than zero only if includeAdhesion is set to true
{
	Real contactsAdhesive=0;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
		if(!I->isReal()) continue;
		JiangPhys* phys = dynamic_cast<JiangPhys*>(I->phys.get());
		if (phys->isAdhesive) {contactsAdhesive += 1;}
	}
	return contactsAdhesive;
}

/***************************************************************************************************************/
/* FUNCTION WHICH RETURNS THE RATIO BETWEEN THE NUMBER OF SLIDING CONTACTS TO THE TOTAL NUMBER AT A GIVEN TIME */
/***************************************************************************************************************/

Real Law2_ScGeom_JiangPhys_Jiang::ratioSlidingContacts()
{
	Real ratio(0); int count(0);
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
		if(!I->isReal()) continue;
		JiangPhys* phys = dynamic_cast<JiangPhys*>(I->phys.get());
		if (phys->isSliding) {ratio+=1;}
		count++;
	}  
	ratio/=count;
	return ratio;
}

/*********************************************************************/
/* FUNCTION TO GET THE NORMAL ELASTIC POTENTIAL ENERGY OF THE SYSTEM */
/*********************************************************************/

Real Law2_ScGeom_JiangPhys_Jiang::normElastEnergy()
{
	Real normEnergy=0;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
		if(!I->isReal()) continue;
		ScGeom* scg = dynamic_cast<ScGeom*>(I->geom.get());
		JiangPhys* phys = dynamic_cast<JiangPhys*>(I->phys.get());
		if (phys) {
			if (includeAdhesion) {normEnergy += pow((scg->penetrationDepth*phys->kn - phys->adhesionForce*scg->penetrationDepth),2)/(2*phys->kn);}
			else {normEnergy += pow(scg->penetrationDepth*phys->kn,2)/(2*phys->kn);} // work done in the normal direction
			}
	}
	return normEnergy;
}

/*****************************************************/
/* FUNCTION TO GET THE ADHESION ENERGY OF THE SYSTEM */
/*****************************************************/

Real Law2_ScGeom_JiangPhys_Jiang::adhesionEnergy()
{
	Real adhesionEnergy=0;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
		if(!I->isReal()) continue;
		ScGeom* scg = dynamic_cast<ScGeom*>(I->geom.get());
		JiangPhys* phys = dynamic_cast<JiangPhys*>(I->phys.get());
		if (phys && includeAdhesion) {
			Real R=scg->radius1*scg->radius2/(scg->radius1+scg->radius2);
			Real gammapi=phys->adhesionForce/(4.*R);
			adhesionEnergy += gammapi*pow(phys->radius,2);} // note that contact radius is calculated if we calculate energy components (for now radius is from Hertz theory)
	}
	return adhesionEnergy;
}

/**********************************************************/
/************** Law2_ScGeom_JiangPhys_Jiang ***************/
/**********************************************************/

CREATE_LOGGER(Law2_ScGeom_JiangPhys_Jiang);

bool Law2_ScGeom_JiangPhys_Jiang::go(shared_ptr<IGeom>& ig, shared_ptr<IPhys>& ip, Interaction* contact){
	const Real& dt = scene->dt; // get time step
	
	Body::id_t id1 = contact->getId1(); // get id body 1
 	Body::id_t id2 = contact->getId2(); // get id body 2

	State* de1 = Body::byId(id1,scene)->state.get();
	State* de2 = Body::byId(id2,scene)->state.get();	

	ScGeom* scg = static_cast<ScGeom*>(ig.get());
	JiangPhys* phys = static_cast<JiangPhys*>(ip.get());	

	const shared_ptr<Body>& b1=Body::byId(id1,scene); 
	const shared_ptr<Body>& b2=Body::byId(id2,scene); 

	bool useViscDamping=(phys->betanDir!=0. || phys->betasDir!=0. || phys->betanIndir!=0.);
	bool DirViscDamp=true;
	if (phys->betanIndir!=0.) {DirViscDamp=false;} // use indirect viscous damping by defining en and or es

	// tangential,normal,rolling and twisting viscous damping coefficients, recomputed from betan,betas at every step
	Real cn=0, cs=0, cr=0, ct=0;

	/****************/
	/* NORMAL FORCE */
	/****************/
	
	Real uN = scg->penetrationDepth; // get overlapping 
	if (uN<0) {
		if (neverErase) {phys->shearForce = phys->normalForce = Vector3r::Zero(); phys->kn=phys->ks=0; return true;}
		else return false;
	}
	/* Jiang formulation 
	Note that the normal stiffness here is a secant value (so as it is cannot be used in the GSTS)
	In the first place we get the normal force and then we store kn to be passed to the GSTS */
	Real Fn = phys->kn*uN; // normal Force (scalar)
	if (includeAdhesion) {
			Fn -= phys->adhesionForce; // include adhesion force to account for the effect of Van der Waals interactions
			phys->isAdhesive = (Fn<0); // set true the bool to count the number of adhesive contacts
			}
	phys->normalForce = Fn*scg->normal; // normal Force (vector)


	if (calcEnergy){
		Real R=scg->radius1*scg->radius2/(scg->radius1+scg->radius2);
		phys->radius=pow((Fn+(includeAdhesion?phys->adhesionForce:0.))*pow(R,3/2.)/phys->kn,1/3.); // attribute not used anywhere, we do not need it
	}

	
	/********************************/
	/* VISCOUS DAMPING COEFFICIENTS */
	/********************************/
	
	// Inclusion of local damping if requested
	// viscous damping is defined for both direct (assigning betan and betas) and indirect (assigning en and es)
	if (useViscDamping && DirViscDamp){// this means that direct method must be used by assigning betan and betas directly
		Real mbar = (!b1->isDynamic() && b2->isDynamic()) ? de2->mass : ((!b2->isDynamic() && b1->isDynamic()) ? de1->mass : (de1->mass*de2->mass / (de1->mass + de2->mass))); // get equivalent mass if both bodies are dynamic, if not set it equal to the one of the dynamic body
		//Real mbar = de1->mass*de2->mass / (de1->mass + de2->mass); // equivalent mass
		Real Cn_crit = 2.*sqrt(mbar*phys->kn); // critical damping coefficient (normal direction)
		Real Cs_crit = 2.*sqrt(mbar*phys->ks); // critical damping coefficient (shear direction)
		// Note: to compare with the analytical solution you provide cn and cs directly (since here we used a different method to define c_crit)
		cn = Cn_crit*phys->betanDir; // coefficient of viscous normal damping 
		cs = Cs_crit*phys->betasDir; // coefficient of viscous shear  damping
		cr = 0.25*pow(phys->R_bar,2)*cn; // Coefficient of viscous rolling damping (equation 13 of Jiang's paper) 
		ct = 0.5*pow(phys->R_bar,2)*cs; // Coefficient of viscous twisting damping (equation 25 of Jiang's paper)
		if(phys->kn<0 || phys->ks<0){ cerr<<"Negative stiffness kn="<<phys->kn<<" ks="<<phys->ks<<" for ##"<<b1->getId()<<"+"<<b2->getId()<<", step "<<scene->iter<<endl; }
	}
	else if (useViscDamping){ // this means that indirect method must be used by assigning en for calculating betanIndir and then cn and cs
		Real mbar = (!b1->isDynamic() && b2->isDynamic()) ? de2->mass : ((!b2->isDynamic() && b1->isDynamic()) ? de1->mass : (de1->mass*de2->mass / (de1->mass + de2->mass))); // get equivalent mass if both bodies are dynamic, if not set it equal to the one of the dynamic body
		cn = phys->betanIndir*2.0*sqrt(mbar*phys->kn); // normal viscous coefficient (equation 2b of Jiang's paper)
		cs = cn; // for now,same value for coefficient of viscous shear damping is assumed. It is subject to further investigation to employ a different formula for cs
		cr = 0.25*pow(phys->R_bar,2)*cn; // Coefficient of viscous rolling damping (equation 13 of Jiang's paper) 
		ct = 0.5*pow(phys->R_bar,2)*cs; // Coefficient of viscous twisting damping (equation 25 of Jiang's paper)
	}


	/***************/
	/* SHEAR FORCE */
	/***************/
	
	Vector3r& shearElastic = phys->shearElastic; // reference for shearElastic force
	// Define shifts to handle periodicity
	const Vector3r shift2 = scene->isPeriodic ? scene->cell->intrShiftPos(contact->cellDist): Vector3r::Zero(); 
	const Vector3r shiftVel = scene->isPeriodic ? scene->cell->intrShiftVel(contact->cellDist): Vector3r::Zero(); 
	// 1. Rotate shear force
	shearElastic = scg->rotate(shearElastic);
	Vector3r prev_FsElastic = shearElastic; // save shear force at previous time step
#if 0 
    // This is my first code for calculating  incidentVs. In this code, I considered rotParticleI and rotParticleJ but it seems that angVelParticleI and angVelParticleJ must be used instead. 
	// 2. Get incident velocity, get shear and normal components
	//see equation 32b of Jiang's paper for figuring out how incident shear velocity is calcualted in our case
	Vector3r incidentV = scg->getIncidentVel(de1, de2, dt, shift2, shiftVel, preventGranularRatcheting);
	Vector3r incidentVn = scg->normal.dot(incidentV)*scg->normal; // contact normal velocity (see equation 32a of Jiang's paper)
	Vector3r angVelParticleI = b1->state->angVel;
	Vector3r angVelParticleJ = b2->state->angVel;
	Vector3r rotParticleI = angVelParticleI*dt;
	Vector3r rotParticleJ = angVelParticleJ*dt;
	Real r_primeI = scg->radius1-(scg->penetrationDepth/2);
	Real r_primeJ = scg->radius2-(scg->penetrationDepth/2);
	Vector3r incidentVs = incidentV - incidentVn + ((r_primeI*scg->normal).cross(rotParticleI)) + ((r_primeJ*scg->normal).cross(rotParticleJ)); // contact shear velocity. This is slightly different from original HM (see equation 32b of Jiang's paper)
#endif	
	// 2. Get incident velocity, get shear and normal components
	//see equation 32b of Jiang's paper for figuring out how incident shear velocity is calcualted in our case
	Vector3r incidentV = scg->getIncidentVel(de1, de2, dt, shift2, shiftVel, preventGranularRatcheting);
	Vector3r incidentVn = scg->normal.dot(incidentV)*scg->normal; // contact normal velocity (see equation 32a of Jiang's paper)
	Vector3r angVelParticleI = b1->state->angVel;
	Vector3r angVelParticleJ = b2->state->angVel;
	Real r_primeI = scg->radius1-(scg->penetrationDepth/2);
	Real r_primeJ = scg->radius2-(scg->penetrationDepth/2);
	Vector3r incidentVs = incidentV - incidentVn + ((r_primeI*scg->normal).cross(angVelParticleI)) + ((r_primeJ*scg->normal).cross(angVelParticleJ)); // contact shear velocity. This is slightly different from original HM (see equation 32b of Jiang's paper)
	// 3. Get shear force (incrementally)
	shearElastic = shearElastic - phys->ks*(incidentVs*dt);



	/***************************************************************************************************************************************/
	/* UPDATING NORMAL FORCE BEFORE MOHR-COLUMB AND MOHR-COLUMB LIKE CRITEIA APPLICATION AND CALCULATING NORMAL VISCOUS DAMPING DISSIPATION*/
	/***************************************************************************************************************************************/
	
	// normal force must be updated here before we apply the Mohr-Coulomb (in shear) and Mohr-Coulumb like (in rolling and twisting) criteria
	if (useViscDamping){ // get normal viscous component
		phys->normalViscous = cn*incidentVn;
		Vector3r normTemp = phys->normalForce - phys->normalViscous; // temporary normal force
		// viscous force should not exceed the value of current normal force, i.e. no attraction force should be permitted if particles are non-adhesive
		// if particles are adhesive, then fixed the viscous force at maximum equal to the adhesion force
		// *** enforce normal force to zero if no adhesion is permitted ***
		if (phys->adhesionForce==0.0 || !includeAdhesion){
						if (normTemp.dot(scg->normal)<0.0){
										phys->normalForce = Vector3r::Zero();
										phys->normalViscous = phys->normalViscous + normTemp; // normal viscous force is such that the total applied force is null - it is necessary to compute energy correctly!
						}
						else{phys->normalForce -= phys->normalViscous;}
		}
		else if (includeAdhesion && phys->adhesionForce!=0.0){
						// *** limit viscous component to the max adhesive force ***
						if (normTemp.dot(scg->normal)<0.0 && (phys->normalViscous.norm() > phys->adhesionForce) ){
										Real normVisc = phys->normalViscous.norm(); Vector3r normViscVector = phys->normalViscous/normVisc;
										phys->normalViscous = phys->adhesionForce*normViscVector;
										phys->normalForce -= phys->normalViscous;
						}
						// *** apply viscous component - in the presence of adhesion ***
						else {phys->normalForce -= phys->normalViscous;}
		}
		if (calcEnergy) {normDampDissip += phys->normalViscous.dot(incidentVn*dt);} // calc dissipation of energy due to normal damping
	}


	/*************************************/
	/* SHEAR DISPLACEMENT (ELASTIC ONLY) */
	/*************************************/
	
	Vector3r& us_elastic = phys->usElastic;
	us_elastic = scg->rotate(us_elastic); // rotate vector
	Vector3r prevUs_el = us_elastic; // store previous elastic shear displacement (already rotated)
	us_elastic -= incidentVs*dt; // add shear increment


	/****************************************/
	/* SHEAR DISPLACEMENT (ELASTIC+PLASTIC) */
	/****************************************/
	
	Vector3r& us_total = phys->usTotal;
	us_total = scg->rotate(us_total); // rotate vector
	Vector3r prevUs_tot = us_total; // store previous total shear displacement (already rotated)
	us_total -= incidentVs*dt; // add shear increment NOTE: this vector is not passed into the failure criterion, hence it holds also the plastic part of the shear displacement

	bool noShearDamp = false; // bool to decide whether we need to account for shear damping dissipation or not


	/**************************************/
	/* MOHR-COULOMB LAW FOR SHAER SLIDING */
	/**************************************/

	phys->isSliding=false;
	phys->shearViscous=Vector3r::Zero(); // reset so that during sliding, the previous values is not there
	Fn = phys->normalForce.norm();
	if (!includeAdhesion) {
		Real maxFs = Fn*phys->tangensOfFrictionAngle;
		if (shearElastic.squaredNorm() > maxFs*maxFs){
			phys->isSliding=true;
			noShearDamp = true; // no damping is added in the shear direction, hence no need to account for shear damping dissipation
			Real ratio = maxFs/shearElastic.norm();
			shearElastic *= ratio; phys->shearForce = shearElastic; /*store only elastic shear displacement*/ us_elastic*= ratio;
			if (calcEnergy) {frictionDissipation += (us_total-prevUs_tot).dot(shearElastic);} // calcualte energy dissipation due to plastic state (sliding) in tangential direction
			}
		else if (useViscDamping){ // add current contact damping if we do not slide and if damping is requested
			phys->shearViscous = cs*incidentVs; // get shear viscous component
			phys->shearForce = shearElastic - phys->shearViscous;}
		else if (!useViscDamping) {phys->shearForce = shearElastic;} // update the shear force at the elastic value if no damping is present and if we passed MC
	}
	else { // Mohr-Coulomb formulation adpated due to the presence of adhesion (see Thornton, 1991).
		Real maxFs = phys->tangensOfFrictionAngle*(phys->adhesionForce+Fn); // adhesionForce already included in normalForce (above)
		if (shearElastic.squaredNorm() > maxFs*maxFs){
			phys->isSliding=true;
			noShearDamp = true; // no damping is added in the shear direction, hence no need to account for shear damping dissipation
			Real ratio = maxFs/shearElastic.norm(); shearElastic *= ratio; phys->shearForce = shearElastic; /*store only elastic shear displacement*/ us_elastic *= ratio;
			if (calcEnergy) {frictionDissipation += (us_total-prevUs_tot).dot(shearElastic);} // calcualte energy dissipation due to plastic state (sliding) in tangential direction
			}
		else if (useViscDamping){ // add current contact damping if we do not slide and if damping is requested
			phys->shearViscous = cs*incidentVs; // get shear viscous component
			phys->shearForce = shearElastic - phys->shearViscous;}
		else if (!useViscDamping) {phys->shearForce = shearElastic;} // update the shear force at the elastic value if no damping is present and if we passed MC
	}


	/**********************************/
	/* SHEAR ELASTIC POTENTIAL ENERGY */
	/**********************************/
	
	// NOTE: shear elastic energy calculation must come after the MC criterion, otherwise displacements and forces are not updated
	if (calcEnergy) {
		shearEnergy += (us_elastic-prevUs_el).dot((shearElastic+prev_FsElastic)/2.); // NOTE: no additional energy if we perform sliding since us_elastic and prevUs_el will hold the same value (in fact us_elastic is only keeping the elastic part). We work out the area of the trapezium.
	}


	/***************************************************/
	/* ENERGY DISSIPATION DUE TO SHEAR VISCOUS DAMPING */
	/***************************************************/
	
	if (useViscDamping){ // get normal viscous component (the shear one is calculated inside Mohr-Coulomb criterion, see above)
		if (calcEnergy) {if (!noShearDamp) {shearDampDissip += phys->shearViscous.dot(incidentVs*dt);}} // calc energy dissipation due to shear viscous damping
	}


	/****************/
	/* APPLY FORCES */
	/****************/
	
	if (!scene->isPeriodic)
		applyForceAtContactPoint(-phys->normalForce - phys->shearForce, scg->contactPoint , id1, de1->se3.position, id2, de2->se3.position);
	else { // in scg we do not wrap particles positions, hence "applyForceAtContactPoint" cannot be used
		Vector3r force = -phys->normalForce - phys->shearForce;
		scene->forces.addForce(id1,force);
		scene->forces.addForce(id2,-force);
		scene->forces.addTorque(id1,(scg->radius1-0.5*scg->penetrationDepth)* scg->normal.cross(force));
		scene->forces.addTorque(id2,(scg->radius2-0.5*scg->penetrationDepth)* scg->normal.cross(force));
	}


	/*****************************/
	/* ROLLING RESISTANCE MOMENT */
	/*****************************/

	if (includeRollResistMoment){
		Vector3r& rollMomentElastic = phys->rollMomentElastic; // reference for rollMomentElastic moment
		// 1. Rotate rolling moment
		rollMomentElastic = scg->rotate(rollMomentElastic); // rotate rolling moment vector (updated)
		Vector3r prev_MrElastic = rollMomentElastic; // save rolling moment at previous time step
		// 2. Compute relative particle rotation in rolling direction (similar to the way the shear is computed)
		// use scg function to compute relAngVel
		Vector3r relAngVel = scg->getRelAngVel(de1,de2,dt);
		//Vector3r relAngVel = (b2->state->angVel-b1->state->angVel); // an alternative way
		Vector3r relAngVelBend = relAngVel - scg->normal.dot(relAngVel)*scg->normal; // keep only the bending part (equation 33a of Jiang's Paper) 
		Vector3r relRot = relAngVelBend*dt; // relative rotation due to rolling behaviour	
		// 3. Implement incremental formulation for the rolling moment (as for the shear part)
		rollMomentElastic = rollMomentElastic-phys->kr*relRot; // add incremental rolling to the rolling vector
#if 0
	// code to compute the relative particle rotation
	if (includeRollResistMoment){
		Real rMean = (scg->radius1+scg->radius2)/2.;
		// sliding motion
		Vector3r duS1 = scg->radius1*(phys->prevNormal-scg->normal);
		Vector3r duS2 = scg->radius2*(scg->normal-phys->prevNormal);
		// rolling motion
		Vector3r duR1 = scg->radius1*dt*b1->state->angVel.cross(scg->normal);
		Vector3r duR2 = -scg->radius2*dt*b2->state->angVel.cross(scg->normal);
		// relative position of the old contact point with respect to the new one
		Vector3r relPosC1 = duS1+duR1;
		Vector3r relPosC2 = duS2+duR2;
		
		Vector3r duR = (relPosC1+relPosC2)/2.; // incremental displacement vector (same radius is temporarily assumed)

		// check wheter rolling will be present, if not do nothing
		Vector3r x=scg->normal.cross(duR);
		Vector3r normdThetaR(Vector3r::Zero()); // initialize 
		if(x.squaredNorm()==0) { /* no rolling */ }
		else {
				Vector3r normdThetaR = x/x.norm(); // moment unit vector
				phys->dThetaR = duR.norm()/rMean*normdThetaR;} // incremental rolling
		
		// incremental formulation for the bending moment (as for the shear part)
		Vector3r& rollMoment = phys->rollMoment;
		rollMomentElastic = scg->rotate(rollMomentElastic); // rotate moment vector
		rollMomentElastic = rollMomentElastic+phys->kr*phys->dThetaR; // add incremental rolling to the rolling vector FIXME: is the sign correct?
#endif
    /**********************************************************************************************************************************************************************/                                                                                                                                                                  
	/* NOTE: NORMAL FORCE HAS ALREADY BEEN UPDATED BEFORE IMPLEMENTING MOHR-COLUMB CRITERION IN SHEAR. NO NEED TO UPDATE AGAIN FOR MC LIKE CRITERION IN ROLLING DIRECTION */                                                                                                                                                                  
    /**********************************************************************************************************************************************************************/                                                                                                                                                                  
                                                                                                                                                                                                                                                                                                                                              
	                                                                                                                                                                                                                                                                                                                                            
	                                                                                                                                                                                                                                                                                                                                            
	/***********************************/                                                                                                                                                                                                                                                                                                     
    /* ROLLING ROTATION (ELASTIC ONLY) */                                                                                                                                                                                                                                                                                                     
    /***********************************/                                                                                                                                                                                                                                                                                                     
                                                                                                                                                                                                                                                                                                                                              
	  Vector3r& tetaRollElastic = phys->tetaRollElastic;                                                                                                                                                                                                                                                                                        
	  tetaRollElastic = scg->rotate(tetaRollElastic); // rotate vector                                                                                                                                                                                                                                                                          
	  Vector3r prevTetaRollElastic = tetaRollElastic; // store previous elastic rolling rotation (already rotated)                                                                                                                                                                                                                              
	  tetaRollElastic -= relAngVelBend*dt; // add rolling rotation increment                                                                                                                                                                                                                                                                    
	                                                                                                                                                                                                                                                                                                                                            
	                                                                                                                                                                                                                                                                                                                                            
	/**************************************/                                                                                                                                                                                                                                                                                                  
	/* ROLLING ROTATION (ELASTIC+PLASTIC) */                                                                                                                                                                                                                                                                                                  
    /**************************************/                                                                                                                                                                                                                                                                                                  
                                                                                                                                                                                                                                                                                                                                              
	  Vector3r& tetaRollTotal = phys->tetaRollTotal;                                                                                                                                                                                                                                                                                            
	  tetaRollTotal = scg->rotate(tetaRollTotal); // rotate vector                                                                                                                                                                                                                                                                              
	  Vector3r prevTetaRollTotal = tetaRollTotal; // store previous total rolling rotation (already rotated)                                                                                                                                                                                                                                    
	  tetaRollTotal -= relAngVelBend*dt; // add rolling rotation increment NOTE: this vector is not passed into the failure criterion, hence it holds also the plastic part of the rolling rotation                                                                                                                                             
	                                                                                                                                                                                                                                                                                                                                            
	  bool noRollDamp = false; // bool to decide whether we need to account for rolling damping dissipation or not                                                                                                                                                                                                                              
	                                                                                                                                                                                                                                                                                                                                            
	                                                                                                                                                                                                                                                                                                                                            
    /*****************************************************************************************/                                                                                                                                                                                                                                               
	/* MOHR-COULOMB LIKE LAW FOR LIMITING ROLLING RESISTANCE MOMENT AFTER THE ELASTIC REGIME */                                                                                                                                                                                                                                               
    /*****************************************************************************************/                                                                                                                                                                                                                                               
                                                                                                                                                                                                                                                                                                                                              
	  phys->therRollRotReached=false;                                                                                                                                                                                                                                                                                                           
	  phys->rollViscous=Vector3r::Zero(); // reset so that after therRollRotReached, the previous values is not there                                                                                                                                                                                                                           
	  Fn = phys->normalForce.norm();                                                                                                                                                                                                                                                                                                            
    	if (!includeAdhesion){                                                                                                                                                                                                                                                                                                                  
	  		Real maxMr = 0.25*phys->XIC*Fn*phys->R_bar; // see equation 19b and figure 7c in Jiang's paper                                                                                                                                                                                                                                        
	  		if (rollMomentElastic.squaredNorm() > maxMr*maxMr){                                                                                                                                                                                                                                                                                   
	  			phys->therRollRotReached=true;                                                                                                                                                                                                                                                                                                      
	  			noRollDamp = true; // no damping is added in the rolling direction, hence no need to account for rolling damping dissipation                                                                                                                                                                                                        
	  			Real ratio = maxMr/rollMomentElastic.norm();                                                                                                                                                                                                                                                                                        
	  			rollMomentElastic *= ratio; phys->rollMoment = rollMomentElastic; /*store only elastic rolling rotation*/ tetaRollElastic*= ratio;                                                                                                                                                                                                  
	  			if (calcEnergy) {rollingPlasticDissipation += (tetaRollTotal-prevTetaRollTotal).dot(rollMomentElastic);} // calcualte energy dissipation due to plastic state in rolling direction                                                                                                                                                  
	  		else if (useViscDamping){ // add current contact damping if therRollRotReached is not reached and if damping is requested	                                                                                                                                                                                                            
	  			phys->rollViscous = cr*relAngVelBend; // get rolling viscous component                                                                                                                                                                                                                                                              
	  			phys->rollMoment = rollMomentElastic - phys->rollViscous;}                                                                                                                                                                                                                                                                          
	  		else if (!useViscDamping) {phys->rollMoment = rollMomentElastic;} // update the rolling moment at the elastic value if no damping is present and if we passed MC-like law                                                                                                                                                             
	  	}                                                                                                                                                                                                                                                                                                                                       
	  	else { // Mohr-Coulomb like formulation adpated due to the presence of adhesion (see Thornton, 1991).	                                                                                                                                                                                                                                  
	  		maxMr = 0.25*phys->XIC*phys->R_bar*(phys->adhesionForce+Fn); // adhesionForce already included in normalForce (above)                                                                                                                                                                                                            
	  		if (rollMomentElastic.squaredNorm() > maxMr*maxMr){                                                                                                                                                                                                                                                                                   
	  			phys->therRollRotReached=true;                                                                                                                                                                                                                                                                                                      
	  			noRollDamp = true; // no damping is added in the rolling direction, hence no need to account for rolling damping dissipation                                                                                                                                                                                                        
	  			Real ratio = maxMr/rollMomentElastic.norm();                                                                                                                                                                                                                                                                                        
	  			rollMomentElastic *= ratio; phys->rollMoment = rollMomentElastic; /*store only elastic rolling rotation*/ tetaRollElastic*= ratio;                                                                                                                                                                                                  
	  			if (calcEnergy) {rollingPlasticDissipation += (tetaRollTotal-prevTetaRollTotal).dot(rollMomentElastic);} // calcualte energy dissipation due to plastic state in rolling direction                                                                                                                                                  
	  		    }                                                                                                                                                                                                                                                                                                                                 
	                                                                                                                                                                                                                                                                                                                                            
	  		else if (useViscDamping){ // add current contact damping if therRollRotReached is not reached and if damping is requested	                                                                                                                                                                                                            
	  			phys->rollViscous = cr*relAngVelBend; // get rolling viscous component                                                                                                                                                                                                                                                              
	  			phys->rollMoment = rollMomentElastic - phys->rollViscous;}                                                                                                                                                                                                                                                                          
    		else if (!useViscDamping) {phys->rollMoment = rollMomentElastic;} // update the rolling moment at the elastic value if no damping is present and if we passed MC-like law                                                                                                                                                             
	  	}                                                                                                                                                                                                                                                                                                                                       
	                                                                                                                                                                                                                                                                                                                                           
	                                                                                                                                                                                                                                                                                                                                            
	                                                                                                                                                                                                                                                                                                                                            
	/************************************/                                                                                                                                                                                                                                                                                                    
	/* ROLLING ELASTIC POTENTIAL ENERGY */                                                                                                                                                                                                                                                                                                    
    /************************************/                                                                                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                                                                                              
	  // NOTE: Rolling elastic energy calculation must come after the MC like criterion, otherwise displacements and forces are not updated                                                                                                                                                                                                     
	  if (calcEnergy) {                                                                                                                                                                                                                                                                                                                         
	  	rollEnergy += (tetaRollElastic-prevTetaRollElastic).dot((rollMomentElastic+prev_MrElastic)/2.); // NOTE: no additional energy if therRollRotReached is reached since tetaRollElastic and prevTetaRollElastic will hold the same value (in fact tetaRollElastic is only keeping the elastic part). We work out the area of the trapezium.
	  }                                                                                                                                                                                                                                                                                                                                         
	                                                                                                                                                                                                                                                                                                                                            
	                                                                                                                                                                                                                                                                                                                                            
	/*****************************************************/                                                                                                                                                                                                                                                                                   
	/* ENERGY DISSIPATION DUE TO ROLLING VISCOUS DAMPING */                                                                                                                                                                                                                                                                                   
    /*****************************************************/                                                                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                                                                                                              
	  if (useViscDamping){                                                                                                                                                                                                                                                                                                                      
	  	if (calcEnergy) {if (!noRollDamp) {rollingDampDissip += phys->rollViscous.dot(relAngVelBend*dt);}} // calc energy dissipation due to rolling viscous damping                                                                                                                                                                            
	  } 
	}                                                                                                                                                                                                                                                                                                                                          
	
	
	/******************************/
	/* TWISTING RESISTANCE MOMENT */
	/******************************/

	if (includeTwistResistMoment){
		Vector3r& twistMomentElastic = phys->twistMomentElastic; // reference for twistMomentElastic moment
		// 1. Rotate twisting moment
		twistMomentElastic = scg->rotate(twistMomentElastic); // rotate twisting moment vector (updated)
		Vector3r prev_MtElastic = twistMomentElastic; // save twisting moment at previous time step
		// 2. Compute relative particle rotation in twisting direction (similar to the way the shear is computed)
		// use scg function to compute relAngVel
		relAngVel = (b2->state->angVel-b1->state->angVel);
		Vector3r relAngVelTwist = scg->normal.dot(relAngVel)*scg->normal; // (equation 33b of Jiang's Paper)
		Vector3r relRotTwist = relAngVelTwist*dt; // component of relative rotation along n
		// 3. Implement incremental formulation for the twisting moment
		twistMomentElastic = twistMomentElastic-phys->ktw*relRotTwist; // add incremental twisting to the twisting vector  
		
	/**********************************************************************************************************************************************************************/                                                                                                                                                                              
	/* NOTE: NORMAL FORCE HAS ALREADY BEEN UPDATED BEFORE IMPLEMENTING MOHR-COLUMB CRITERION IN SHEAR. NO NEED TO UPDATE AGAIN FOR MC LIKE CRITERION IN ROLLING DIRECTION */                                                                                                                                                                              
    /**********************************************************************************************************************************************************************/                                                                                                                                                                              
                                                                                                                                                                                                                                                                                                                                                          
	                                                                                                                                                                                                                                                                                                                                                        
	/************************************/                                                                                                                                                                                                                                                                                                                
	/* TWISTING ROTATION (ELASTIC ONLY) */                                                                                                                                                                                                                                                                                                                
    /************************************/                                                                                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                                                                                                                          
	  Vector3r& tetaTwistElastic = phys->tetaTwistElastic;                                                                                                                                                                                                                                                                                                  
	  tetaTwistElastic = scg->rotate(tetaTwistElastic); // rotate vector                                                                                                                                                                                                                                                                                    
	  Vector3r prevTetaTwistElastic = tetaTwistElastic; // store previous elastic twisting rotation (already rotated)                                                                                                                                                                                                                                       
	  tetaTwistElastic -= relAngVelTwist*dt; // add twisting rotation increment                                                                                                                                                                                                                                                                             
	                                                                                                                                                                                                                                                                                                                                                        
	                                                                                                                                                                                                                                                                                                                                                        
	/***************************************/                                                                                                                                                                                                                                                                                                             
	/* TWISTING ROTATION (ELASTIC+PLASTIC) */                                                                                                                                                                                                                                                                                                             
    /***************************************/                                                                                                                                                                                                                                                                                                             
                                                                                                                                                                                                                                                                                                                                                          
	  Vector3r& tetaTwistTotal = phys->tetaTwistTotal;                                                                                                                                                                                                                                                                                                      
	  tetaTwistTotal = scg->rotate(tetaTwistTotal); // rotate vector                                                                                                                                                                                                                                                                                        
	  Vector3r prevTetaTwistTotal = tetaTwistTotal; // store previous total twisting rotation (already rotated)                                                                                                                                                                                                                                             
	  tetaTwistTotal -= relAngVelTwist*dt; // add twisting rotation increment NOTE: this vector is not passed into the failure criterion, hence it holds also the plastic part of the twisting rotation                                                                                                                                                     
	                                                                                                                                                                                                                                                                                                                                                        
	  bool noTwistDamp = false; // bool to decide whether we need to account for twisting damping dissipation or not	                                                                                                                                                                                                                                      
	                                                                                                                                                                                                                                                                                                                                                        
	                                                                                                                                                                                                                                                                                                                                                        
    /******************************************************************************************/                                                                                                                                                                                                                                                          
	/* MOHR-COULOMB LIKE LAW FOR LIMITING TWISTING RESISTANCE MOMENT AFTER THE ELASTIC REGIME */                                                                                                                                                                                                                                                          
    /******************************************************************************************/                                                                                                                                                                                                                                                          
                                                                                                                                                                                                                                                                                                                                                          
	  phys->therTwistRotReached=false;                                                                                                                                                                                                                                                                                                                      
	  phys->twistViscous=Vector3r::Zero(); // reset so that after therTwistRotReached , the previous values is not there                                                                                                                                                                                                                                    
	  Fn = phys->normalForce.norm();                                                                                                                                                                                                                                                                                                                        
    	if (!includeAdhesion){                                                                                                                                                                                                                                                                                                                              
	  		Real maxMt = 0.65*phys->tangensOfFrictionAngle*Fn*phys->R_bar; // see equation 31 and figure 7d in Jiang's paper                                                                                                                                                                                                                                  
	  		if (twistMomentElastic.squaredNorm() > maxMt*maxMt){                                                                                                                                                                                                                                                                                              
	  			phys->therTwistRotReached=true;                                                                                                                                                                                                                                                                                                                 
	  			noTwistDamp = true; // no damping is added in the twisting direction, hence no need to account for twisting damping dissipation                                                                                                                                                                                                                 
	  			Real ratio = maxMt/twistMomentElastic.norm();                                                                                                                                                                                                                                                                                                   
	  			twistMomentElastic *= ratio; phys->twistMoment = twistMomentElastic; /*store only elastic twisting rotation*/ tetaTwistElastic*= ratio;                                                                                                                                                                                                         
	  			if (calcEnergy) {twistPlasticDissipation += (tetaTwistTotal-prevTetaTwistTotal).dot(twistMomentElastic);} // calculate energy dissipation due to plastic state in twisting direction                                                                                                                                                            
	  			}                                                                                                                                                                                                                                                                                                                                               
	  		else if (useViscDamping){ // add current contact damping if therTwistRotReached is not reached and if damping is requested	                                                                                                                                                                                                                      
	  			phys->twistViscous = ct*relAngVelTwist; // get twisting viscous component                                                                                                                                                                                                                                                                       
	  			phys->twistMoment = twistMomentElastic - phys->twistViscous;}                                                                                                                                                                                                                                                                                   
	  		else if (!useViscDamping) {phys->twistMoment = twistMomentElastic;} // update the twisting moment at the elastic value if no damping is present and if we passed MC-like law                                                                                                                                                                      
	  	}                                                                                                                                                                                                                                                                                                                                                   
	  	else { // Mohr-Coulomb like formulation adpated due to the presence of adhesion (see Thornton, 1991).	                                                                                                                                                                                                                                              
	  		Real maxMt = 0.65*phys->tangensOfFrictionAngle*phys->R_bar*(phys->adhesionForce+Fn); // adhesionForce already included in normalForce (above)                                                                                                                                                                                                     
	  		if (twistMomentElastic.squaredNorm() > maxMt*maxMt){                                                                                                                                                                                                                                                                                              
	  			phys->therTwistRotReached=true;                                                                                                                                                                                                                                                                                                                 
	  			noTwistDamp = true; // no damping is added in the twisting direction, hence no need to account for twisting damping dissipation                                                                                                                                                                                                                 
	  			Real ratio = maxMt/twistMomentElastic.norm();                                                                                                                                                                                                                                                                                                   
	  			twistMomentElastic *= ratio; phys->twistMoment = twistMomentElastic; /*store only elastic twisting rotation*/ tetaTwistElastic*= ratio;                                                                                                                                                                                                         
	  			if (calcEnergy) {twistPlasticDissipation += (tetaTwistTotal-prevTetaRollTotal).dot(twistMomentElastic);} // calcualte energy dissipation due to plastic state in twisting direction                                                                                                                                                             
	  			}                                                                                                                                                                                                                                                                                                                                               
	  		else if (useViscDamping){ // add current contact damping if therTwistRotReached is not reached and if damping is requested	                                                                                                                                                                                                                      
	  			phys->twistViscous = ct*relAngVelTwist; // get twisting viscous component                                                                                                                                                                                                                                                                       
	  			phys->twistMoment = twistMomentElastic - phys->twistViscous;}                                                                                                                                                                                                                                                                                   
	  		else if (!useViscDamping) {phys->twistMoment = twistMomentElastic;} // update the twisting moment at the elastic value if no damping is present and if we passed MC-like law                                                                                                                                                                      
	  	}                                                                                                                                                                                                                                                                                                                                                   
	                                                                                                                                                                                                                                                                                                                                                       
	                                                                                                                                                                                                                                                                                                                                                        
	                                                                                                                                                                                                                                                                                                                                                        
	/*************************************/                                                                                                                                                                                                                                                                                                               
	/* TWISTING ELASTIC POTENTIAL ENERGY */                                                                                                                                                                                                                                                                                                               
    /*************************************/                                                                                                                                                                                                                                                                                                               
                                                                                                                                                                                                                                                                                                                                                          
	  // NOTE: tWISTING elastic energy calculation must come after the MC like criterion, otherwise displacements and forces are not updated                                                                                                                                                                                                                
	  if (calcEnergy) {                                                                                                                                                                                                                                                                                                                                     
	  	twistEnergy += (tetaTwistElastic-prevTetaTwistElastic).dot((twistMomentElastic+prev_MtElastic)/2.); // NOTE: no additional energy if therTwistRotReached is reached since tetaTwistElastic and prevtetaTwistElastic will hold the same value (in fact tetaTwistElastic is only keeping the elastic part). We work out the area of the trapezium.    
	  }                                                                                                                                                                                                                                                                                                                                                     
	                                                                                                                                                                                                                                                                                                                                                        
	                                                                                                                                                                                                                                                                                                                                                        
	/******************************************************/                                                                                                                                                                                                                                                                                              
	/* ENERGY DISSIPATION DUE TO TWISTING VISCOUS DAMPING */                                                                                                                                                                                                                                                                                              
    /******************************************************/                                                                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                                                                                                                                          
	  if (useViscDamping){                                                                                                                                                                                                                                                                                                                                  
	  	if (calcEnergy) {if (!noTwistDamp) {twistDampDissip += phys->twistViscous.dot(relAngVelTwist*dt);}} // calc energy dissipation due to twisting viscous damping                                                                                                                                                                                      
	  }
	}
	
	
	/*****************/
	/* APPLY MOMENTS */
	/*****************/

		Vector3r moment = phys->twistMoment+phys->rollMoment;
		scene->forces.addTorque(id1,-moment); 
		scene->forces.addTorque(id2,moment);
}
return true;
}


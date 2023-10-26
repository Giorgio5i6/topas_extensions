// Scorer for AmorphousTrack
//
// ********************************************************************
// *                                                                  *
// *                                                                  *
// * This file was obtained from Topas MC Inc under the license       *
// * agreement set forth at http://www.topasmc.org/registration       *
// * Any use of this file constitutes full acceptance of              *
// * this TOPAS MC license agreement.                                 *
// *                                                                  *
// ********************************************************************
//

#include "TsTrackerHit_AT.hh"
#include "MyScoreAmorphousTrack.hh"
#include<vector>
#include<fstream>
#include <unordered_map>
#include<string>

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4EmCalculator.hh"
#include "G4UIcommand.hh"

MyScoreAmorphousTrack::MyScoreAmorphousTrack(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
		G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
	: TsVBinnedScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)
{

	G4String lutfolder = fPm->GetStringParameter(GetFullParmName("LUTFolder"));
	foutfilename = fPm->GetStringParameter(GetFullParmName("OutputFile"));
	
	fTypeOfTable = "zdD"; 
	if ( fPm->ParameterExists(GetFullParmName("TypeOfLUT")) ) {
		G4int index = fPm->GetIntegerParameter(GetFullParmName("TypeOfLUT"));
		if(index == 1) fTypeOfTable = "zdDstar";
		if(index == 2) fTypeOfTable = "znD";
	}



	//PARAMETERS
	zStart = 0.1 ;  // in unit of keV/um
	if ( fPm->ParameterExists(GetFullParmName("SpecificEnergyLowerLimit")) ) {
		zStart  =  fPm->GetUnitlessParameter(GetFullParmName("SpecificEnergyLowerLimit"));
	}

	zEnd = 100 ;  // in unit of keV/um
	if ( fPm->ParameterExists(GetFullParmName("SpecificEnergyUpperLimit")) ) {
		zEnd  =  fPm->GetUnitlessParameter(GetFullParmName("SpecificEnergyUpperLimit"));
	}

	fCountsLowerLimit = 100;
	if ( fPm->ParameterExists(GetFullParmName("HistoCountsLowerLimit")) ) {
		fCountsLowerLimit  =  fPm->GetIntegerParameter(GetFullParmName("HistoCountsLowerLimit"));
	}


	SetUnit("");	
	fHitsCollection = new TsTrackerHitsCollection();

	flut_zD.clear();
	GetLUT(flut_zD, lutfolder, fTypeOfTable);


	fzdDSpectraBin.clear();
	fzdDDoseAverzdD.clear(); // Accumulated for all steps and events

	//INITIALIZE LOG BIN HISTO
	zBins = 100;
	//zStart = 0.1;
	//zEnd = 100;
	zBinLimit.resize(zBins+1, 0.);
	zBinWidth.resize(zBins, 0.);
	zBinCenter.resize(zBins, 0.);

	zBinLimit[0] = zStart; //lowest z value
	double zMax = zEnd; //highest z value
	double step = (log10(zMax) - log10(zStart))/zBins;

	double Binlog = log10(zStart);
	for (int i=0; i<zBins; i++)
	{
		Binlog += step;
		zBinLimit[i+1] = pow(10, Binlog);
		zBinWidth[i] = zBinLimit[i+1]-zBinLimit[i];
		zBinCenter[i] = (zBinLimit[i+1]+zBinLimit[i])/2.;
	}

	fElectronDefinition = G4ParticleTable::GetParticleTable()->FindParticle("e-");
	fStepCount = 0;
}


MyScoreAmorphousTrack::~MyScoreAmorphousTrack() {;}


G4bool MyScoreAmorphousTrack::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
	if (!fIsActive) {
		fSkippedWhileInactive++;
		return false;
	}


	G4double edep = aStep->GetTotalEnergyDeposit();
	if (edep==0.) return false;

	TsTrackerHit* newHit = new TsTrackerHit();

	//G4double weight = aStep->GetPreStepPoint()->GetWeight();
	G4double ekin =  aStep->GetPreStepPoint()->GetKineticEnergy();
	//G4String ParticleName =aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName();
	G4int Z = aStep->GetTrack()->GetParticleDefinition()->GetAtomicNumber();
	G4int A = aStep->GetTrack()->GetParticleDefinition()->GetAtomicMass();
	if(Z>8) return false;

	//Include endep from secondary electrons
	const G4TrackVector* secondary = aStep->GetSecondary();
	if (!secondary) {
		secondary = aStep->GetTrack()->GetStep()->GetSecondary();  // parallel worlds
	}
	if (secondary) {
		G4int diff;
		if (fStepCount == 0) diff = 0;
		else diff = secondary->size() - fStepCount;

		fStepCount = (*secondary).size();

		for (unsigned int i=(*secondary).size()-diff; i<(*secondary).size(); i++)
			if ((*secondary)[i]->GetParticleDefinition() == fElectronDefinition)
				edep += (*secondary)[i]->GetKineticEnergy();
	}


	newHit->SetEkin(ekin);
	newHit->SetBinIndex(GetIndex(aStep));
	newHit->SetParticleAtmMass(A);
	newHit->SetParticleCharge(Z);
	newHit->SetEdep(edep);
	fHitsCollection->insert(newHit);
	return true;
}

void MyScoreAmorphousTrack::AccumulateEvent()
{


	std::unordered_map<G4int, std::vector<G4double>> bin_edepzD;
	G4int nofHits = fHitsCollection->entries();
	for ( G4int i=0; i<nofHits; i++ )
	{

		G4double edep = (*fHitsCollection)[i]->GetEdep();
		G4double ekin = (*fHitsCollection)[i]->GetEkin();
		G4int index = (*fHitsCollection)[i]->GetBinIndex();
		G4int Z =  (*fHitsCollection)[i]->GetParticleCharge();
		G4int A =  (*fHitsCollection)[i]->GetParticleAtmMass();

		//CALCULTAE HISTO f(z_domain)
		std::vector<G4double> zdD(3,0.);
		if(Z == 1) zdD = Interpolate(flut_zD.at("H"),  ekin/A, 1);
		else if(Z == 2) zdD = Interpolate(flut_zD.at("He"), ekin/A, 1);
		//else if(Z == 3) zdD = Interpolate(flut_zD.at("Li"), ekin/A, 1);
		//else if(Z == 4) zdD = Interpolate(flut_zD.at("Be"), ekin/A, 1);
		//else if(Z == 5) zdD = Interpolate(flut_zD.at("B"),ekin/A, 1);
		//else if(Z == 6) zdD = Interpolate(flut_zD.at("C"),ekin/A, 1);
		//else if(Z == 7) zdD = Interpolate(flut_zD.at("N"),ekin/A, 1);
		//else if(Z == 8) zdD = Interpolate(flut_zD.at("O"),ekin/A, 1);


		//Dynamically append at the end of the array if the hit belongs to a new voxel
		if(bin_edepzD.find(index) == bin_edepzD.end())
			bin_edepzD.insert({{index, {edep, edep*zdD[0], edep*zdD[1], edep*zdD[2]}}});
		//If the voxel is already registered
		else
		{
			(bin_edepzD.at(index))[0] += edep; //denominator of the ratio sum(edep*zD)/sum(edep)
			for(int ztype=0; ztype<3; ztype++)
				(bin_edepzD.at(index))[ztype+1] += edep*zdD[ztype]; //numerator of the ratio sum(edep*zD)/sum(edep)
		}

	}



	//Event-by-Event Histo with dose-averaged zdD
	//AT THE MOMENT NOT USED COULD BE REMOVED TO SPEED UP
	for (auto& x: bin_edepzD)
	{

		//x.first is the ID of the voxel
		if(fzdDSpectraBin.find(x.first) == fzdDSpectraBin.end())
		{	
			std::vector<G4double> histo(zBins,0.);
			FillHisto( ((x.second)[1])/((x.second)[0]), histo);
			fzdDSpectraBin.insert({{x.first, histo}});
			fzdDDoseAverzdD.insert({{x.first, x.second}});
		}
		else
		{
			FillHisto( ((x.second)[1])/((x.second)[0]), fzdDSpectraBin.at(x.first));
			for(int ztype=0; ztype<4; ztype++)
				(fzdDDoseAverzdD.at(x.first))[ztype] += (x.second)[ztype];

		}
	}

	delete fHitsCollection;
	fHitsCollection = new TsTrackerHitsCollection();
}

void MyScoreAmorphousTrack::AbsorbResultsFromWorkerScorer(TsVScorer* workerScorer) {
	TsVBinnedScorer::AbsorbResultsFromWorkerScorer(workerScorer);


	MyScoreAmorphousTrack* workerMTScorer = dynamic_cast<MyScoreAmorphousTrack*>(workerScorer);
	std::unordered_map<int, std::vector<double>> fzdDSpectraBin_worker = workerMTScorer->fzdDSpectraBin;
	std::unordered_map<int, std::vector<double>> fzdDDoseAverzdD_worker = workerMTScorer->fzdDDoseAverzdD;

	for (auto& x: fzdDSpectraBin_worker)
	{
		if(fzdDSpectraBin.find(x.first) == fzdDSpectraBin.end())
		{	
			fzdDSpectraBin.insert({{x.first, x.second}});
			fzdDDoseAverzdD.insert({x.first, fzdDDoseAverzdD_worker.at(x.first)});
		}
		else
		{
			for (int i=0; i<zBins; i++)
				(fzdDSpectraBin.at(x.first))[i] += x.second[i];

			for(int ztype=0; ztype<4; ztype++)
				(fzdDDoseAverzdD.at(x.first))[ztype] += (fzdDDoseAverzdD_worker.at(x.first))[ztype];
		}
	}

	workerMTScorer->fzdDSpectraBin.clear();
	workerMTScorer->fzdDDoseAverzdD.clear();

}

void MyScoreAmorphousTrack::UserHookForEndOfRun()
{

	/*
	   std::ofstream outfile("prova.csv");
	   for (auto& x: fzdDSpectraBin)
	   {
	   int sum = 0;
	   for(auto bin:x.second)
	   sum += bin;
	   if(sum > 1)
	   {
	   outfile << x.first;
	   G4cout << x.first;
	   for(auto bin:x.second)
	   {
	   outfile <<','<<bin;
	   G4cout <<','<<bin;
	   }
	   outfile << std::endl;
	   G4cout << std::endl;
	   }
	   }

	   outfile.close();
	 */
}

void MyScoreAmorphousTrack::FillHisto(G4double data, std::vector<G4double> &hist)
{

	for (int n=0;n<zBins;n++)
	{
		if(data<=zBinLimit[n+1])
		{
			hist[n] = hist[n]+1;
			break;
		}
	}
}



std::vector<G4double> MyScoreAmorphousTrack::Interpolate(std::vector<std::vector<G4double>> Data, G4double x, G4bool extrapolate )
{
	std::vector<G4double> xData = Data[0];
	std::vector<std::vector<G4double>> yData = {Data[1], Data[2], Data[3]};
	G4int size = xData.size();

	G4int i = 0;                                                                  // find left end of interval for interpolation
	if ( x >= xData[size - 2] )                                                 // special case: beyond right end
	{
		i = size - 2;
	}
	else
	{
		while ( x > xData[i+1] ) i++;
	}
	G4double xL = xData[i], xR= xData[i+1];
	G4double yL0 = yData[0][i], yR0 = yData[0][i+1];      // points on either side (unless beyond ends)
	G4double yL1 = yData[1][i], yR1 = yData[1][i+1];      // points on either side (unless beyond ends)
	G4double yL2 = yData[2][i], yR2 = yData[2][i+1];      // points on either side (unless beyond ends)
	if ( !extrapolate )                                                         // if beyond ends of array and not extrapolating
	{
		if ( x < xL ) 
		{
			yR0 = yL0;
			yR1 = yL1;
			yR2 = yL2;
		}
		if ( x > xR ) 
		{
			yL0 = yR0;
			yL1 = yR1;
			yL2 = yR2;
		}
	}

	G4double dydx0 = ( yR0 - yL0 ) / ( xR - xL );                                    // gradient
	G4double dydx1 = ( yR1 - yL1 ) / ( xR - xL );                                    // gradient
	G4double dydx2 = ( yR2 - yL2 ) / ( xR - xL );                                    // gradient

	return {yL0 + dydx0 * ( x - xL ), yL1 + dydx1 * ( x - xL ), yL2 + dydx2 * ( x - xL )};                                              // linear interpolation
}

void MyScoreAmorphousTrack::GetLUT(std::unordered_map<std::string, std::vector<std::vector<G4double>> > & lut, std::string lutfolder, std::string ztype)
{

	//std::vector<std::string> ListOfIons = {"H", "He", "Li", "Be", "B", "C", "N", "O"};
	std::vector<std::string> ListOfIons = {"H", "He"};
	for(auto ionType:ListOfIons)
	{
		std::string filename = lutfolder + "/LUT_"+ ionType + ".dat";
		std::ifstream lutfile(filename);

		std::vector<double> lut_ekin, lut_zdD, lut_zdS, lut_znD;

		if(lutfile.fail()) // checks to see if file opended 
		{
			std::cout << "ERROR::FILE NOT FOUND!!!\t" << filename << std::endl;
		}

		G4cout << "\nREADING FILE: " << filename << G4endl << "Ekin/A \t z\n";
		while(!lutfile.eof()) // reads file to end of *file*, not line
		{ 

			double x, y0, y1, y2;
			lutfile >> x 
				>> y0 >> y1 >> y2;


			lut_ekin.push_back(x);
			lut_zdD.push_back(y0);
			lut_zdS.push_back(y1);
			lut_znD.push_back(y2);
			G4cout << x <<'\t' << y0 <<'\t' <<y1 <<'\t'<<y2 << G4endl;
		}
		lutfile.close();

		lut.insert({ionType, {lut_ekin, lut_zdD, lut_zdS, lut_znD}});
	}
}


void MyScoreAmorphousTrack::Output()
{

	G4int run = GetRunID();
	G4String runFileName = foutfilename+'-'+std::to_string(run)+".csv";
	std::ofstream outfile(runFileName.c_str());

	G4String runFileNameSpectra = foutfilename+'-'+std::to_string(run)+"_zdD_Spectra.csv";
	std::ofstream outfileSpectra(runFileNameSpectra.c_str());
	for (auto& idx:fzdDSpectraBin)
	{
		int TotalCounts = 0;
		for(auto zSpectraCounts:idx.second)
			TotalCounts += zSpectraCounts;

		if(TotalCounts > fCountsLowerLimit)
		{
			outfile << idx.first;
			for(int ztype=1; ztype<4; ztype++)
			{
				double DoseAveragezd =  (fzdDDoseAverzdD.at(idx.first))[ztype]/ (fzdDDoseAverzdD.at(idx.first))[0];
				outfile << ',' <<DoseAveragezd;
			}

			outfile << std::endl;
			
			outfileSpectra << idx.first;
			for(auto zSpectraCounts:idx.second)
			{
				outfileSpectra <<','<<zSpectraCounts;
			}
			outfileSpectra << std::endl;
		}
	}

	outfile.close();
	outfileSpectra.close();
}

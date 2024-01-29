#ifndef CODE_SCRIPT_TEMPLATE_H
#define CODE_SCRIPT_TEMPLATE_H

#include "script.h"

namespace Code {

class ScriptTemplate {
public:
  constexpr ScriptTemplate(const char * name, const char * value) : m_name(name), m_value(value) {}
  // INSERT FROM OUTPUT.TXT HERE

static const ScriptTemplate * ActiveAndPassiveImmunization();

static const ScriptTemplate * AdaptiveResponse();

static const ScriptTemplate * AntibodyFunction();

static const ScriptTemplate * AntigenicVariation();

static const ScriptTemplate * AntigenPresentingLeukocytes();

static const ScriptTemplate * AntigenRecognitionBCells();

static const ScriptTemplate * AntigenRecognitionTCells();

static const ScriptTemplate * AntimicrobalPeptidesAndProteins();

static const ScriptTemplate * BarrierDefenses();

static const ScriptTemplate * BCellsAfterSelection();

static const ScriptTemplate * BCellsAntibodies();

static const ScriptTemplate * BTCellsDevelopment();

static const ScriptTemplate * BTCellsDiversity();

static const ScriptTemplate * CellMediatedResponse();

static const ScriptTemplate * CellularAerobicRespiration();

static const ScriptTemplate * CellularInnateDefense();

static const ScriptTemplate * Compartmentalization();

static const ScriptTemplate * ConsequencesOfGeneMutations();

static const ScriptTemplate * CytotoxicTCells();

static const ScriptTemplate * Definitions();

static const ScriptTemplate * DifferentialGeneExpression();

static const ScriptTemplate * DisruptionsInImmuneSystem();

static const ScriptTemplate * DNAElectrophoresis();

static const ScriptTemplate * DNAToRNA();

static const ScriptTemplate * Epidemiology();

static const ScriptTemplate * ExportSynthesizedProtein();

static const ScriptTemplate * Flu();

static const ScriptTemplate * GeneMutations();

static const ScriptTemplate * GeneRegulation();

static const ScriptTemplate * GenomeMutations();

static const ScriptTemplate * HelperTCells();

static const ScriptTemplate * Heritability();

static const ScriptTemplate * HIV();

static const ScriptTemplate * HumoralResponse();

static const ScriptTemplate * ImmuneRejection();

static const ScriptTemplate * ImmunologicalMemory();

static const ScriptTemplate * InflammatoryResponses();

static const ScriptTemplate * InnateImmunityOfInvertebrates();

static const ScriptTemplate * InnateImmunityOfVertebrates();

static const ScriptTemplate * Interleukins();

static const ScriptTemplate * LactoseOperon();

static const ScriptTemplate * LymphaticSystem();

static const ScriptTemplate * MechanismusOfPostTranscriptionalRegulation();

static const ScriptTemplate * MemoryBCells();

static const ScriptTemplate * MemoryTCells();

static const ScriptTemplate * MHCI();

static const ScriptTemplate * MHCII();

static const ScriptTemplate * NoncodingRNAInControllingGeneExpression();

static const ScriptTemplate * Operons();

static const ScriptTemplate * OxidizedReduced();

static const ScriptTemplate * Penetrance();

static const ScriptTemplate * Phagocytosis();

static const ScriptTemplate * Photosynthesis();

static const ScriptTemplate * PlasmaCells();

static const ScriptTemplate * PositiveGeneRegulation();

static const ScriptTemplate * ProliferationBTCells();

static const ScriptTemplate * ProtonMotiveForce();

static const ScriptTemplate * PublicHealth();

static const ScriptTemplate * RecognitionAndResponse();

static const ScriptTemplate * RegulationOfTranscriptionInitiation();

static const ScriptTemplate * SelfToleranceBTCells();

static const ScriptTemplate * TCellsAfterSelection();

static const ScriptTemplate * Transcription1st();

static const ScriptTemplate * Translation2nd();

static const ScriptTemplate * TryptophanOperon();

static const ScriptTemplate * TypesOfChromosomeMutations();

static const ScriptTemplate * TypesOfMutations();
  // ---------------------------------------------------------------------------------------------

  static const ScriptTemplate * Empty();
  const char * name() const { return m_name; }
  const char * content() const { return m_value + Script::StatusSize(); }
  const char * value() const { return m_value; }
private:
  const char * m_name;
  const char * m_value; // holds the 'importation status' and 'current importation status' flags concatenated with the script content
};

}

#endif

#include "script_store.h"

namespace Code {

constexpr char ScriptStore::k_scriptExtension[];

bool ScriptStore::ScriptNameIsFree(const char * baseName) {
  return ScriptBaseNamed(baseName).isNull();
}

ScriptStore::ScriptStore() {
  // INSERT FROM OUTPUT.TXT HERE

addScriptFromTemplate(ScriptTemplate::ActiveAndPassiveImmunization());

addScriptFromTemplate(ScriptTemplate::AdaptiveResponse());

addScriptFromTemplate(ScriptTemplate::AntibodyFunction());

addScriptFromTemplate(ScriptTemplate::AntigenicVariation());

addScriptFromTemplate(ScriptTemplate::AntigenPresentingLeukocytes());

addScriptFromTemplate(ScriptTemplate::AntigenRecognitionBCells());

addScriptFromTemplate(ScriptTemplate::AntigenRecognitionTCells());

addScriptFromTemplate(ScriptTemplate::AntimicrobalPeptidesAndProteins());

addScriptFromTemplate(ScriptTemplate::BarrierDefenses());

addScriptFromTemplate(ScriptTemplate::BCellsAfterSelection());

addScriptFromTemplate(ScriptTemplate::BCellsAntibodies());

addScriptFromTemplate(ScriptTemplate::BTCellsDevelopment());

addScriptFromTemplate(ScriptTemplate::BTCellsDiversity());

addScriptFromTemplate(ScriptTemplate::CellMediatedResponse());

addScriptFromTemplate(ScriptTemplate::CellularAerobicRespiration());

addScriptFromTemplate(ScriptTemplate::CellularInnateDefense());

addScriptFromTemplate(ScriptTemplate::Compartmentalization());

addScriptFromTemplate(ScriptTemplate::ConsequencesOfGeneMutations());

addScriptFromTemplate(ScriptTemplate::CytotoxicTCells());

addScriptFromTemplate(ScriptTemplate::Definitions());

addScriptFromTemplate(ScriptTemplate::DifferentialGeneExpression());

addScriptFromTemplate(ScriptTemplate::DisruptionsInImmuneSystem());

addScriptFromTemplate(ScriptTemplate::DNAElectrophoresis());

addScriptFromTemplate(ScriptTemplate::DNAToRNA());

addScriptFromTemplate(ScriptTemplate::Epidemiology());

addScriptFromTemplate(ScriptTemplate::ExportSynthesizedProtein());

addScriptFromTemplate(ScriptTemplate::Flu());

addScriptFromTemplate(ScriptTemplate::GeneMutations());

addScriptFromTemplate(ScriptTemplate::GeneRegulation());

addScriptFromTemplate(ScriptTemplate::GenomeMutations());

addScriptFromTemplate(ScriptTemplate::HelperTCells());

addScriptFromTemplate(ScriptTemplate::Heritability());

addScriptFromTemplate(ScriptTemplate::HIV());

addScriptFromTemplate(ScriptTemplate::HumoralResponse());

addScriptFromTemplate(ScriptTemplate::ImmuneRejection());

addScriptFromTemplate(ScriptTemplate::ImmunologicalMemory());

addScriptFromTemplate(ScriptTemplate::InflammatoryResponses());

addScriptFromTemplate(ScriptTemplate::InnateImmunityOfInvertebrates());

addScriptFromTemplate(ScriptTemplate::InnateImmunityOfVertebrates());

addScriptFromTemplate(ScriptTemplate::Interleukins());

addScriptFromTemplate(ScriptTemplate::LactoseOperon());

addScriptFromTemplate(ScriptTemplate::LymphaticSystem());

addScriptFromTemplate(ScriptTemplate::MechanismusOfPostTranscriptionalRegulation());

addScriptFromTemplate(ScriptTemplate::MemoryBCells());

addScriptFromTemplate(ScriptTemplate::MemoryTCells());

addScriptFromTemplate(ScriptTemplate::MHCI());

addScriptFromTemplate(ScriptTemplate::MHCII());

addScriptFromTemplate(ScriptTemplate::NoncodingRNAInControllingGeneExpression());

addScriptFromTemplate(ScriptTemplate::Operons());

addScriptFromTemplate(ScriptTemplate::OxidizedReduced());

addScriptFromTemplate(ScriptTemplate::Penetrance());

addScriptFromTemplate(ScriptTemplate::Phagocytosis());

addScriptFromTemplate(ScriptTemplate::Photosynthesis());

addScriptFromTemplate(ScriptTemplate::PlasmaCells());

addScriptFromTemplate(ScriptTemplate::PositiveGeneRegulation());

addScriptFromTemplate(ScriptTemplate::ProliferationBTCells());

addScriptFromTemplate(ScriptTemplate::ProtonMotiveForce());

addScriptFromTemplate(ScriptTemplate::PublicHealth());

addScriptFromTemplate(ScriptTemplate::RecognitionAndResponse());

addScriptFromTemplate(ScriptTemplate::RegulationOfTranscriptionInitiation());

addScriptFromTemplate(ScriptTemplate::SelfToleranceBTCells());

addScriptFromTemplate(ScriptTemplate::TCellsAfterSelection());

addScriptFromTemplate(ScriptTemplate::Transcription1st());

addScriptFromTemplate(ScriptTemplate::Translation2nd());

addScriptFromTemplate(ScriptTemplate::TryptophanOperon());

addScriptFromTemplate(ScriptTemplate::TypesOfChromosomeMutations());

addScriptFromTemplate(ScriptTemplate::TypesOfMutations());


  // -----------------------------------------------------------------
}

void ScriptStore::deleteAllScripts() {
  for (int i = numberOfScripts() - 1; i >= 0; i--) {
    scriptAtIndex(i).destroy();
  }
}

bool ScriptStore::isFull() {
  return Ion::Storage::sharedStorage()->availableSize() < k_fullFreeSpaceSizeLimit;
}

const char * ScriptStore::contentOfScript(const char * name, bool markAsFetched) {
  Script script = ScriptNamed(name);
  if (script.isNull()) {
    return nullptr;
  }
  if (markAsFetched) {
    script.setFetchedFromConsole(true);
  }
  return script.content();
}

void ScriptStore::clearVariableBoxFetchInformation() {
  // TODO optimize fetches
  const int scriptsCount = numberOfScripts();
  for (int i = 0; i < scriptsCount; i++) {
    scriptAtIndex(i).setFetchedForVariableBox(false);
  }
}

void ScriptStore::clearConsoleFetchInformation() {
  // TODO optimize fetches
  const int scriptsCount = numberOfScripts();
  for (int i = 0; i < scriptsCount; i++) {
    scriptAtIndex(i).setFetchedFromConsole(false);
  }
}

Script::ErrorStatus ScriptStore::addScriptFromTemplate(const ScriptTemplate * scriptTemplate) {
  size_t valueSize = Script::StatusSize() + strlen(scriptTemplate->content()) + 1; // (auto importation status + content fetched status) + scriptcontent size + null-terminating char
  assert(Script::nameCompliant(scriptTemplate->name()));
  Script::ErrorStatus err = Ion::Storage::sharedStorage()->createRecordWithFullName(scriptTemplate->name(), scriptTemplate->value(), valueSize);
  assert(err != Script::ErrorStatus::NonCompliantName);
  return err;
}

}

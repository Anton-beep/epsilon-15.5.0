#include "script_store.h"

namespace Code {

constexpr char ScriptStore::k_scriptExtension[];

bool ScriptStore::ScriptNameIsFree(const char * baseName) {
  return ScriptBaseNamed(baseName).isNull();
}

ScriptStore::ScriptStore() {
  // INSERT FROM OUTPUT.TXT HERE

  addScriptFromTemplate(ScriptTemplate::searcher());

addScriptFromTemplate(ScriptTemplate::importer());

addScriptFromTemplate(ScriptTemplate::AlternatingCurrent());

addScriptFromTemplate(ScriptTemplate::BFieldOfCoilOfWireSolenoid());

addScriptFromTemplate(ScriptTemplate::Diffraction());

addScriptFromTemplate(ScriptTemplate::DopplerEffect());

addScriptFromTemplate(ScriptTemplate::ElectricField());

addScriptFromTemplate(ScriptTemplate::ElectricPotential());

addScriptFromTemplate(ScriptTemplate::ElectromagneticInduction());

addScriptFromTemplate(ScriptTemplate::ElectromagneticSpectrum());

addScriptFromTemplate(ScriptTemplate::Electronvolt());

addScriptFromTemplate(ScriptTemplate::FaradayLaw());

addScriptFromTemplate(ScriptTemplate::Generator());

addScriptFromTemplate(ScriptTemplate::GravitationalAcceleration());

addScriptFromTemplate(ScriptTemplate::HallEffect());

addScriptFromTemplate(ScriptTemplate::HygensPrinciple());

addScriptFromTemplate(ScriptTemplate::InterferenceWaves());

addScriptFromTemplate(ScriptTemplate::LenzLaw());

addScriptFromTemplate(ScriptTemplate::LightWaves());

addScriptFromTemplate(ScriptTemplate::LorentzForce());

addScriptFromTemplate(ScriptTemplate::LrOrRlCircuits());

addScriptFromTemplate(ScriptTemplate::MagneticFields());

addScriptFromTemplate(ScriptTemplate::MagneticFlux());

addScriptFromTemplate(ScriptTemplate::MassSpectrometry());

addScriptFromTemplate(ScriptTemplate::MotionalEmf());

addScriptFromTemplate(ScriptTemplate::NumberOfMaxima());

addScriptFromTemplate(ScriptTemplate::OerstedsExperiment());

addScriptFromTemplate(ScriptTemplate::RefractionOfLight());

addScriptFromTemplate(ScriptTemplate::RlcCircuit());

addScriptFromTemplate(ScriptTemplate::SelfInduction());

addScriptFromTemplate(ScriptTemplate::SlitsWaves());

addScriptFromTemplate(ScriptTemplate::SoundWaves());

addScriptFromTemplate(ScriptTemplate::StandingWaves());

addScriptFromTemplate(ScriptTemplate::Superposition());

addScriptFromTemplate(ScriptTemplate::Torque());

addScriptFromTemplate(ScriptTemplate::TypesOfWaves());

addScriptFromTemplate(ScriptTemplate::WaveEquation());

addScriptFromTemplate(ScriptTemplate::WavesFeatures());

addScriptFromTemplate(ScriptTemplate::WavesInRippleTank());

addScriptFromTemplate(ScriptTemplate::WavesGrating());

addScriptFromTemplate(ScriptTemplate::Youngsdoubleslit());

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

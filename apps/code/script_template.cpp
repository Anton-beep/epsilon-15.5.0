#include "script_template.h"

namespace Code {

// INSERT FROM OUTPUT.TXT HERE

constexpr ScriptTemplate ActiveAndPassiveImmunizationScriptTemplate("ActiveAndPassiveImmunization.py", "\x01" R"(Active immunity develops
naturally when memory cells
form clones in response to
an infection. It can also
develop following
immunization, also called
vaccination. In
immunization, a
nonpathogenic form of a
microbe or part of a microbe
elicits an immune response
to an immunological memory.
Passive immunity provides
immediate, short-term
protection. It is conferred
naturally when IgG
(antibody) crosses the
placenta from mother to
fetus or when IgA passes
from mother to infant in
breast milk. It can be
conferred artificially by
injecting antibodies into a
nonimmune person.)");

const ScriptTemplate * ScriptTemplate::ActiveAndPassiveImmunization() {return &ActiveAndPassiveImmunizationScriptTemplate;}

constexpr ScriptTemplate AdaptiveResponseScriptTemplate("AdaptiveResponse.py", "\x01" R"(The adaptive response relies
on two types of lymphocytes,
or white blood cells.
Lymphocytes that mature in
the thymus above the heart
are called T cells and those
that mature in bone marrow
are called B cells.
Antigens are substances that
can elicit a response from a
B or T cell. Exposure to the
pathogen activates B and T
cells with antigen
receptors. The small
accessible part of an
antigen that binds to an
antigen receptor is called
an epitope.
B cells have Y (looks like
Y) on the surface, T cells
have I 
B cells and T cells have
receptor proteins that can
bind to foreign molecules.
Each individual lymphocyte
is specialized to recognize
a specific type of molecule.
see AntigenRecognition... (B
and T cells)

Acquired (adaptive) immunity
has two branches: the
humoral immune response and
the cell mediated immune
response. In the humoral
immune response antibodies
help neutralize or eliminate
toxins and pathogens in the
blood and lymph. In the cell
mediated immune response
specialized T cells destroy
affected host cells.

Both the humoral and cell
mediated responses
can include primary and
secondary immune
response. Memory cells
enable the secondary
response.)");

const ScriptTemplate * ScriptTemplate::AdaptiveResponse() {return &AdaptiveResponseScriptTemplate;}

constexpr ScriptTemplate AntibodyFunctionScriptTemplate("AntibodyFunction.py", "\x01" R"(Antibodies do not kill
pathogens; instead they mark
pathogens for destruction.
In neutralization,
antibodies bind to viral
surface proteins preventing
infection of a host cell. 
Antibodies may also bind to
toxins in body fluids and
prevent them from entering
body cells.
In opsonization, antibodies
bind to antigens on bacteria
creating a target for
macrophages or neutrophils,
triggering phagocytosis.
Antigen antibody complexes
may bind to a complement
protein which triggers a
cascade of complement
protein activation.
Ultimately a membrane attack
complex forms a pore in the
membrane of the foreign
cell, leading to its lysis.

Binding of antibodies to
antigens inactivates
antigens by:
neutralization(blocks viral
binding sites; coats
bacteria), agglutination of
microbes, precipitation of
dissolved antigens. These
enhance phagocytosis by
macrophage.
Binding ... by activation of
complement system leads to
cell lysis.

Antibody specificity and
antigen antibody binding has
been harnessed in research,
diagnosis, and therapy.
Polyclonal antibodies,
produced following exposure
to a microbial antigen, are
products of many different
clones of plasma cells, each
specific for a different
epitope. Monoclonal
antibodies are prepared from
a single clone of B cells
grown in culture.)");

const ScriptTemplate * ScriptTemplate::AntibodyFunction() {return &AntibodyFunctionScriptTemplate;}

constexpr ScriptTemplate AntigenicVariationScriptTemplate("AntigenicVariation.py", "\x01" R"(Through antigenic variation,
some pathogens are able to
change epitope expression
and prevent recognition. 
The human influenza virus
mutates rapidly, and new flu
vaccines must be made each
year. Human viruses
occasionally exchange genes
with the viruses of
domesticated animals. This
poses a danger as human
immune systems are unable to
recognize the new viral
strain. Some viruses may
remain in a host in an
inactive state called
latency.)");

const ScriptTemplate * ScriptTemplate::AntigenicVariation() {return &AntigenicVariationScriptTemplate;}

constexpr ScriptTemplate AntigenPresentingLeukocytesScriptTemplate("AntigenPresentingLeukocytes.py", "\x01" R"(Antigen-presenting
leukocytes, such as
dendritic cells,
macrophages, and B cells,
play a crucial role in
linking the nonspecific
(innate) and adaptive immune
responses. This process
involves two key aspects:
recognition of pathogens or
antigens and activation of
specific immune responses.
First, antigen-presenting
leukocytes recognize and
engulf pathogens or antigens
through various pattern
recognition receptors
(PRRs), such as Toll-like
receptors (TLRs) and
scavenger receptors. This
recognition triggers
nonspecific innate immune
responses, such as
inflammation and activation
of other immune cells, to
limit the spread of
infection and initiate the
clearance of pathogens.
Next, antigen-presenting
leukocytes process and
present fragments of the
engulfed antigens, known as
antigen presentation, on
their cell surfaces using
major histocompatibility
complex (MHC) molecules.
This process allows them to
display antigens to other
immune cells, particularly T
lymphocytes (T cells), which
are essential for adaptive
immune responses.
T cells survey the
antigen-presenting
leukocytes and recognize
specific antigens presented
on MHC molecules through
their T cell receptors
(TCRs). This recognition
leads to the activation and
differentiation of T cells
into effector cells, such as
cytotoxic T cells or helper
T cells, depending on the
context of the antigen
presentation.
Effector T cells then carry
out specific functions to
eliminate pathogens or
modulate immune responses.
For example, cytotoxic T
cells directly kill infected
or abnormal cells, while
helper T cells release
cytokines to activate other
immune cells and enhance the
overall immune response.)");

const ScriptTemplate * ScriptTemplate::AntigenPresentingLeukocytes() {return &AntigenPresentingLeukocytesScriptTemplate;}

constexpr ScriptTemplate AntigenRecognitionBCellsScriptTemplate("AntigenRecognitionBCells.py", "\x01" R"(Each B cell antigen receptor
is a Y shaped molecule with
two identical heavy chains
and two identical light
chains. The constant regions
of the chains vary little
among B cells, whereas the
variable regions differ
greatly. The variable
regions provide antigen
specificity.
Binding of a B cell antigen
receptor to an pathogen's
epitope is an early step in
B cell activation. This
gives rise to cells that
secrete a soluble form of
the protein called an
antibody or immunoglobulin
(Ig). Secreted antibodies
are similar to B cell
receptors but lack
transmembrane regions that
anchor receptors in the
plasma membrane.

Bcells looks like:
Have antigen receptors (look
like Y) on the surface. On
the top there are
Antigen-binding sites,
between broken lines there
is a Disulfide bridge.)");

const ScriptTemplate * ScriptTemplate::AntigenRecognitionBCells() {return &AntigenRecognitionBCellsScriptTemplate;}

constexpr ScriptTemplate AntigenRecognitionTCellsScriptTemplate("AntigenRecognitionTCells.py", "\x01" R"(Each T cell receptor
consists of two different
polypeptide chains (called
alpha and beta). The tips of
the chain form a variable
(V) region; the rest is a
constant (C) region.
T cells bind to antigen
fragments displayed or
presented on a host cell.
These antigen fragments are
bound to cell surface
proteins called MHC
molecules. MHC (major
histocompatibility complex)
molecules are host proteins
that display the antigen
fragments on the cell
surface. In infected cells,
MHC molecules bind and
transport antigen fragments
to the cell surface, a
process called antigen
presentation. A T cell can
then bind both the antigen
fragment and the MHC
molecule. This interaction
is necessary for the T cell
to participate in the
adaptive immune response.

Tcells look like:
Have antigen receptors (look
like I, but in zoom looks
like II) on the surface.
Antigen binding cite on the
top of a receptor. One of II
is alpha chain, second is
beta chain. Between II
disulfide bridge.)");

const ScriptTemplate * ScriptTemplate::AntigenRecognitionTCells() {return &AntigenRecognitionTCellsScriptTemplate;}

constexpr ScriptTemplate AntimicrobalPeptidesAndProteinsScriptTemplate("AntimicrobalPeptidesAndProteins.py", "\x01" R"(Peptides and proteins
function in innate defense
by attacking pathogens or
impeding their reproduction.
Interferon proteins provide
innate defense, interfering
with viruses and helping
activate macrophages. About
30 proteins make up the
complement system, which
causes lysis of invading
cells and helps trigger
inflammation.)");

const ScriptTemplate * ScriptTemplate::AntimicrobalPeptidesAndProteins() {return &AntimicrobalPeptidesAndProteinsScriptTemplate;}

constexpr ScriptTemplate BarrierDefensesScriptTemplate("BarrierDefenses.py", "\x01" R"(Vertebrates have
several BARRIER
DEFENSES that protect them
from pathogens. These
defenses are part of the
body’s most basic defense
mechanisms and are
continuously working to
protect against a broad
range of pathogens. Here are
some examples of barrier
defenses in vertebrates:
- The skin is the primary
barrier to the entrance of
microorganisms into the
body. The skin is covered
with a layer of dead,
keratinized epithelium that
is too dry for bacteria to
grow. Additionally, sweat
and other skin secretions
may lower pH, contain toxic
lipids, and physically wash
microbes away.
- Mucous membranes line the
respiratory, urinary, and
reproductive tracts, and
trap both microbes and
debris, facilitating their
removal. Mucus also contains
various antimicrobial
substances, such as
lysozyme, lactoferrin, and
defensins, which can kill or
inhibit the growth of
microbes.
- Many body fluids, such as
saliva, mucus, and tears,
are hostile to microbes. For
example, the saliva in the
mouth is rich in lysozyme,
an enzyme that destroys
bacteria by digesting their
cell walls. The low pH of
the digestive system also
prevents the growth of
microbes.
- Normal flora refers to the
nonpathogenic bacteria that
colonize the skin, mucous
membranes, and other
surfaces of the body. Normal
flora can compete with
pathogenic bacteria for
nutrients and space, and
produce antimicrobial
substances that inhibit the
growth of pathogens.
- Tears contain lysozyme, an
enzyme that can destroy
bacteria by breaking down
their cell walls. Tears also
help to wash away debris and
irritants from the eyes.
- Earwax is a sticky
substance that traps dust,
dirt, and microbes,
preventing them from
entering the ear
canal. Earwax also contains
antimicrobial substances
that can kill or inhibit the
growth of bacteria and
fungi.
- Urine is a waste product
of the urinary system that
can flush out microbes from
the urethra and bladder. The
low pH of urine also
inhibits the growth of
bacteria.)");

const ScriptTemplate * ScriptTemplate::BarrierDefenses() {return &BarrierDefensesScriptTemplate;}

constexpr ScriptTemplate BCellsAfterSelectionScriptTemplate("BCellsAfterSelection.py", "\x01" R"(Roles: B-cells are primarily
responsible for mediating
humoral immunity by
producing antibodies, also
known as immunoglobulins.
These antibodies can
neutralize pathogens, mark
them for destruction by
other immune cells, or
facilitate their clearance
from the body. B-cells also
play a role in antigen
presentation to helper
T-cells and the regulation
of immune responses.
Origin: B-cells originate
from hematopoietic stem
cells (HSCs) in the bone
marrow. During their
development, B-cell
precursors undergo several
stages of differentiation
and maturation within the
bone marrow
microenvironment, ultimately
giving rise to mature,
antigen-inexperienced
B-cells that migrate to
secondary lymphoid organs,
such as the spleen and lymph
nodes, where they await
encounter with antigens.)");

const ScriptTemplate * ScriptTemplate::BCellsAfterSelection() {return &BCellsAfterSelectionScriptTemplate;}

constexpr ScriptTemplate BCellsAntibodiesScriptTemplate("BCellsAntibodies.py", "\x01" R"(The humoral response is
characterized by secretion
of antibodies by B cells.
Activation of the humoral
immune response involves B
cells and helper T cells as
well as proteins on the
surface of pathogens. In
response to cytokines from
helper T cells and an
antigen, a B cell
proliferates and
differentiates into memory B
cells and antibody secreting
effector cells called plasma
cells.)");

const ScriptTemplate * ScriptTemplate::BCellsAntibodies() {return &BCellsAntibodiesScriptTemplate;}

constexpr ScriptTemplate BTCellsDevelopmentScriptTemplate("BTCellsDevelopment.py", "\x01" R"(After selection, both
B-cells and T-cells
circulate through secondary
lymphoid organs, such as
lymph nodes and the spleen,
where they encounter
antigens presented by
antigen-presenting cells
(APCs), such as dendritic
cells.
Upon encountering specific
antigens, B-cells and
T-cells undergo activation,
proliferation, and
differentiation into
effector cells. B-cells
differentiate into plasma
cells, which produce and
secrete antibodies tailored
to the encountered antigen.
T-cells differentiate into
various effector subsets,
such as helper T-cells or
cytotoxic T-cells, depending
on the nature of the antigen
and the cytokine milieu.
Effector B-cells and T-cells
participate in the immune
response by eliminating
pathogens, coordinating
immune responses, and
promoting the clearance of
infected or abnormal cells.
Following the resolution of
the immune response, a
subset of activated B-cells
and T-cells differentiate
into long-lived memory
cells, which persist in the
body and provide rapid and
enhanced immune responses
upon re-exposure to the same
antigen.)");

const ScriptTemplate * ScriptTemplate::BTCellsDevelopment() {return &BTCellsDevelopmentScriptTemplate;}

constexpr ScriptTemplate BTCellsDiversityScriptTemplate("BTCellsDiversity.py", "\x01" R"(By combining variable
elements, the immune system
assembles a diverse variety
of antigen
receptors. The
immunoglobulin (Ig) gene
encodes one chain of the B
cell receptor. Many
different chains can be
produced from the same gene
by rearrangement of the DNA.
Rearranged DNA is
transcribed and translated
and the antigen receptor
formed.

Somatic Recombination (B
cells): B cells undergo
somatic recombination during
their development in the
bone marrow. This process
involves the rearrangement
of gene segments encoding
the variable regions of the
B cell receptor (BCR),
namely the V (variable), D
(diversity), and J (joining)
gene segments. This random
recombination results in a
vast repertoire of BCRs with
different antigen-binding
specificities.

B cells can express five
different forms (or classes)
of immunoglobulin (Ig) with
similar antigen binding
specificity but different
heavy chain C regions:
- IgD: Membrane bound
- IgM: First soluble class
produced
- IgG: Second soluble class;
most abundant
- IgA and IgE: Remaining
soluble classes)");

const ScriptTemplate * ScriptTemplate::BTCellsDiversity() {return &BTCellsDiversityScriptTemplate;}

constexpr ScriptTemplate CellMediatedResponseScriptTemplate("CellMediatedResponse.py", "\x01" R"(Cell-mediated immune
responses are a vital part
of the body's defense
against intracellular
pathogens, including
viruses, some bacteria,
fungi, and parasites, as
well as cancer cells. Unlike
humoral immunity, which
involves the production of
antibodies by B cells,
cell-mediated immunity
relies on the actions of
specialized white blood
cells called T lymphocytes
or T cells.
During a cell-mediated
immune response,
antigen-presenting cells,
such as dendritic cells,
macrophages, or B cells,
present antigens from
pathogens to T cells. This
presentation activates
certain subsets of T cells,
including cytotoxic T cells
(also known as CD8+ T cells)
and helper T cells (also
known as CD4+ T cells).
Cytotoxic T cells are
primarily responsible for
directly attacking and
killing infected or abnormal
cells. They recognize
specific antigens displayed
on the surface of infected
cells or cancer cells, which
are usually derived from
intracellular pathogens or
abnormal cellular proteins.
Once activated, cytotoxic T
cells release cytotoxic
molecules, such as perforin
and granzymes, which induce
apoptosis (programmed cell
death) in the target cells,
effectively eliminating
them.
Helper T cells play a
crucial role in coordinating
and regulating the immune
response. They release
signaling molecules called
cytokines, which activate
and enhance the functions of
other immune cells,
including cytotoxic T cells,
B cells, and macrophages.
Helper T cells also help
activate memory T cells,
which can provide
long-lasting immunity upon
re-exposure to the same
antigen.
Cell-mediated immune
responses are essential for
clearing intracellular
pathogens and controlling
infections that reside
within host cells. They are
particularly important in
combating viral infections,
as viruses often replicate
inside host cells, evading
detection by antibodies.
Additionally, cell-mediated
immunity plays a critical
role in immune surveillance
against cancer cells,
helping to identify and
eliminate abnormal cells
before they develop into
tumors.
Understanding cell-mediated
immune responses is crucial
for developing vaccines,
immunotherapies, and
treatments for infectious
diseases, cancer, autoimmune
disorders, and other
conditions where T cell
responses are dysregulated.
Manipulating these responses
can help modulate the immune
system to promote protective
immunity or suppress harmful
immune reactions.)");

const ScriptTemplate * ScriptTemplate::CellMediatedResponse() {return &CellMediatedResponseScriptTemplate;}

constexpr ScriptTemplate CellularAerobicRespirationScriptTemplate("CellularAerobicRespiration.py", "\x01" R"(Formula: C6H12O6 + 6O2 ->
6CO2 + 6H2O + 36-38ATP
#### Glycolysis
Glycolysis is a sequence of
chemical reactions that
breaks down a glucose into
two pyruvate molecules. 2
Pyruvate molecules produced
during glycolysis and then
enter the mitochondria.
#### Krebs cycle
The cycle is made up of
eight steps catalyzed by
eight different enzymes that
produce energy at several
different stages.
#### Oxidative
phosphorylation (electron
transport chain)
Electrons are passed from
one member of the transport
chain to another in a series
of redox reactions. Energy
released in these reactions
is captured as a proton
gradient, which is then used
to make ATP in a process
called CHEMIOSMOSIS.
Together, the electron
transport chain and
chemiosmosis make
up oxidative
phosphorylation.
As the electrons travel
through the chain, they go
from a higher to a lower
energy level. Energy is
released in these “downhill”
electron transfers, and
several of the protein
complexes use the released
energy to pump protons from
the mitochondrial matrix to
the intermembrane space,
forming a proton gradient.
Like many other ions,
protons can't pass directly
through the phospholipid
bilayer of the membrane
because its core is too
hydrophobic. Instead,
H+ ions can move down their
concentration gradient only
with the help of channel
proteins that form
hydrophilic tunnels across
the membrane.
In the inner mitochondrial
membrane, H+ ions have just
one channel available: a
membrane-spanning protein
known as ATP (ADENOSINE
TRIPHOSPHATE) SYNTHASE.
Conceptually, ATP synthase
is a lot like a turbine in a
hydroelectric power plant.
Instead of being turned by
water, it’s turned by the
flow of H+ ions moving down
their electrochemical
gradient. As ATP synthase
turns, it catalyzes the
addition of a phosphate to
ADP, capturing energy from
the proton gradient as ATP.)");

const ScriptTemplate * ScriptTemplate::CellularAerobicRespiration() {return &CellularAerobicRespirationScriptTemplate;}

constexpr ScriptTemplate CellularInnateDefenseScriptTemplate("CellularInnateDefense.py", "\x01" R"(Phagocytic cells recognize
groups of pathogens by TLRs
(Toll-like receptors)
A white blood cell engulfs a
microbe, then fuses with a
lysosome to destroy the
microbe.
Types of phagocytic cells:
- Neutrophils: engulf and
destroy pathogens
- Macrophages: are found
throughout the body
- Dendritic cells stimulate
development of adaptive
immunity
- Eosinophils discharge
destructive enzymes
Cellular innate defenses in
vertebrates also involve
natural killer cells. These
circulate through the body
and detect abnormal cells.
They release chemicals
leading to cell death,
inhibiting the spread of
virally infected or
cancerous cells. Many
cellular innate defenses
involve the lymphatic
system.)");

const ScriptTemplate * ScriptTemplate::CellularInnateDefense() {return &CellularInnateDefenseScriptTemplate;}

constexpr ScriptTemplate CompartmentalizationScriptTemplate("Compartmentalization.py", "\x01" R"(Compartmentalization refers
to the way organelles in
eukaryotic cells live and
work in separate areas
within the cell in order to
perform their specific
functions more efficiently.
Eukaryotic cells contain a
membrane-bound nucleus, as
well as membrane-bound
organelles that each perform
different functions. Those
organelles live within
different compartments
inside the cell, so they can
work in the microenvironment
that suits them best or use
separations to create a
gradient between areas.
In mitochondria:
intermembrane space/matrix,
in the chloroplast:
thylakoid/ lumen/stroma).)");

const ScriptTemplate * ScriptTemplate::Compartmentalization() {return &CompartmentalizationScriptTemplate;}

constexpr ScriptTemplate ConsequencesOfGeneMutationsScriptTemplate("ConsequencesOfGeneMutations.py", "\x01" R"(- mRNA: Gene mutations can
alter the sequence of mRNA
transcribed from the DNA,
which can affect the codons
that specify the amino acids
in the protein. Depending on
the type and location of the
mutation, this can result in
different outcomes, such as
same-sense mutations (no
change in amino acid),
missense mutations (change
in amino acid), nonsense
mutations (premature stop
codon), or frameshift
mutations (insertion or
deletion of nucleotides that
shifts the reading frame).
- Ribosomes: Gene mutations
also affect the structure
and function of ribosomes,
which are the molecular
machines that translate mRNA
into proteins. Ribosomes are
composed of ribosomal RNA
(rRNA) and proteins, and
mutations in either
component can impair the
assembly, stability, or
activity of the ribosome.
For example, mutations in
rRNA genes can cause defects
in ribosome biogenesis,
ribosomal accuracy, or
antibiotic resistance.
- tRNA: Gene mutations can
influence the processing and
function of transfer RNA
(tRNA), which are the
molecules that carry amino
acids to the ribosome and
pair with the mRNA codons.
Mutations in tRNA genes or
tRNA processing enzymes can
affect the charging,
folding, or recognition of
tRNA, leading to reduced
translation efficiency,
mistranslation, or diseases.
For instance, mutations in
mitochondrial tRNA (mt-tRNA)
genes are associated with
various human disorders,
such as diabetes, deafness,
or epilepsy.
- Proteins: Gene mutations
can ultimately affect the
structure and function of
proteins, which are the
products of gene expression
and perform various
biological roles. Mutations
in protein-coding genes can
alter the amino acid
sequence, folding,
stability, interactions, or
activity of proteins,
resulting in
loss-of-function,
gain-of-function, or
dominant-negative effects.
Depending on the protein
involved, mutations can
cause phenotypic changes,
diseases, or evolutionary
adaptations.)");

const ScriptTemplate * ScriptTemplate::ConsequencesOfGeneMutations() {return &ConsequencesOfGeneMutationsScriptTemplate;}

constexpr ScriptTemplate CytotoxicTCellsScriptTemplate("CytotoxicTCells.py", "\x01" R"(Cytotoxic T cells are
responsible for identifying
and destroying infected or
abnormal cells.
Functions:
1. Recognition of Infected
Cells: Cytotoxic T cells are
activated when they
encounter cells displaying
foreign antigens on their
surface. These antigens are
typically derived from
pathogens that have invaded
the host cell or from
abnormal proteins expressed
by cancer cells. The
recognition of these
antigens is facilitated by T
cell receptors (TCRs) on the
surface of cytotoxic T
cells, which bind
specifically to antigenic
peptides presented in the
context of major
histocompatibility complex
(MHC) class I molecules on
the surface of target cells.
2. Activation and Expansion:
Upon activation, cytotoxic T
cells undergo clonal
expansion, proliferating to
generate a large population
of effector cells capable of
recognizing and eliminating
infected or abnormal cells.
This activation process is
initiated by signals from
antigen-presenting cells
(APCs), such as dendritic
cells, which present
antigenic peptides to
cytotoxic T cells and
provide co-stimulatory
signals necessary for their
activation.
3. Effector Functions:
Activated cytotoxic T cells
exert their effector
functions through various
mechanisms to eliminate
target cells. One of the
primary mechanisms is the
release of cytotoxic
molecules, such as perforin
and granzymes, which induce
apoptosis (programmed cell
death) in the target cell.
Perforin forms pores in the
target cell's membrane,
allowing granzymes to enter
and trigger apoptotic
pathways. Additionally,
cytotoxic T cells can
express cell surface
proteins, such as Fas ligand
(FasL), which can engage
with Fas receptors on target
cells, leading to apoptosis.
4. Memory and Longevity:
Like other T cell subsets,
cytotoxic T cells can
differentiate into memory
cells following the
resolution of an infection.
Memory cytotoxic T cells
persist in the body for
extended periods, providing
long-term immunity against
recurrent infections by
mounting rapid and robust
responses upon re-exposure
to the same antigen.
5. Tumor Surveillance: In
addition to their role in
combating infectious
pathogens, cytotoxic T cells
are involved in immune
surveillance against cancer.
They can recognize and
eliminate cancerous cells
that express tumor-specific
antigens, helping to prevent
the development and spread
of tumors.)");

const ScriptTemplate * ScriptTemplate::CytotoxicTCells() {return &CytotoxicTCellsScriptTemplate;}

constexpr ScriptTemplate DefinitionsScriptTemplate("Definitions.py", "\x01" R"(DROSOMYCIN is an antifungal
peptide from Drosophila
melanogaster (fruit fly)
DEFENSIN is an antibacterial
peptide with activity
against Gram positive
bacteria.)");

const ScriptTemplate * ScriptTemplate::Definitions() {return &DefinitionsScriptTemplate;}

constexpr ScriptTemplate DifferentialGeneExpressionScriptTemplate("DifferentialGeneExpression.py", "\x01" R"(Almost all the cells in an
organism are genetically
identical. Differences
between cell types result
from differential gene
expression, the expression
of different genes by cells
with the same genome.
Abnormalities in gene
expression can lead to
diseases including cancer.
Gene expression is regulated
at many stages.

Histone Modifications: In
histone acetylation, acetyl
groups are attached to
positively charged lysines
in histone tails. This
loosens chromatin structure,
thereby promoting the
initiation of transcription.
The addition of methyl
groups (methylation) can
condense chromatin; the
addition of phosphate groups
(phosphorylation) next to a
methylated amino acid can
loosen chromatin.

DNA Methylation: DNA
methylation, the addition of
methyl groups to certain
bases in DNA, is associated
with reduced transcription
in some species. DNA
methylation can cause
long-term inactivation of
genes in cellular
differentiation. In genomic
imprinting, methylation
regulates expression of
either the maternal or
paternal alleles of certain
genes at the start of
development.

Epigenetic Inheritance:
Although the chromatin
modifications just discussed
do not alter DNA sequence,
they may be passed to future
generations of cells. The
inheritance of traits
transmitted by mechanisms
not directly involving the
nucleotide sequence is
called epigenetic
inheritance.)");

const ScriptTemplate * ScriptTemplate::DifferentialGeneExpression() {return &DifferentialGeneExpressionScriptTemplate;}

constexpr ScriptTemplate DisruptionsInImmuneSystemScriptTemplate("DisruptionsInImmuneSystem.py", "\x01" R"(Allergies:
Allergies are
exaggerated (hypersensitive)
responses to antigens called
allergens. In localized
allergies such as hay fever,
IgE antibodies produced
after first exposure to an
allergen attach to receptors
on mast cells. The next time
the allergen enters the
body, it binds to mast cell
associated IgE molecules.
Mast cells release histamine
and other mediators that
cause vascular changes
leading to typical allergy
symptoms. An acute allergic
response can lead to
anaphylactic shock, a life
threatening reaction, within
seconds of allergen
exposure.

Autoimmune disease:
In individuals with
autoimmune diseases, the
immune system loses
tolerance for self and turns
against certain molecules of
the body. Autoimmune
diseases include systemic
lupus erythematosus,
rheumatoid arthritis,
insulin dependent diabetes
mellitus, and multiple
sclerosis.

Immunodeficiency disease:
Inborn immunodeficiency
results from hereditary or
developmental defects that
prevent proper functioning
of innate, humoral, and/or
cell mediated defenses.
Acquired immunodeficiency
develops later in life and
results from exposure to
chemical and biological
agents. Acquired
immunodeficiency syndrome
(AIDS) is caused by a virus.)");

const ScriptTemplate * ScriptTemplate::DisruptionsInImmuneSystem() {return &DisruptionsInImmuneSystemScriptTemplate;}

constexpr ScriptTemplate DNAElectrophoresisScriptTemplate("DNAElectrophoresis.py", "\x01" R"(It is a technique to
separate DNA fragments based
on their length. DNA samples
are loaded on a gel and the
gel is connected to a power
supply. Since DNA is
negatively charged due to
the phosphate groups, the
DNA fragments will migrate
to the + pole. Small
fragments will move more
easily through the gel and
end up further from the well
than the larger fragments.
Normal beta-globin cell: 175
- 201 - large
Sickle-cell mutant
beta-globin allele: 376 -
large
HbS - Sickle-cell hemoglobin
HbA - A hemoglobin)");

const ScriptTemplate * ScriptTemplate::DNAElectrophoresis() {return &DNAElectrophoresisScriptTemplate;}

constexpr ScriptTemplate DNAToRNAScriptTemplate("DNAToRNA.py", "\x01" R"(G -> C
C -> G
A -> U
T -> A)");

const ScriptTemplate * ScriptTemplate::DNAToRNA() {return &DNAToRNAScriptTemplate;}

constexpr ScriptTemplate EpidemiologyScriptTemplate("Epidemiology.py", "\x01" R"(It is the study of the
distribution and
determinants of
health-related states and
events in specified
populations, and the
application of this study to
the control of health
problems.
Objectives: Identify risk
factors, Explain the
etiology/cause of disease,
Determine the prognoses of
the disease, Investigate the
effectiveness of
interventions and health
measures (prevention and
treatment).

Endemic is when a disease is
always present in a
particular community or
region. (West Nile virus,
flu, Lyme in NY)
Epidemic is an outbreak of a
disease that affects many
people in one community or
region. (Yellow fever,
smallpox). An epidemic
doesn't have to be caused by
a contagious disease. For
example, heart disease,
obesity and opioid deaths
are also considered
epidemics.
Pandemic is n outbreak of a
disease that has spread to
multiple countries or
continents. (COVID-19)

Systemic diseases: Systemic
means affecting the entire
body, rather than a single
organ or body part. An
infection that is in the
bloodstream is called a
systemic infection. E.g. of
systemic diseases: High
blood pressure, flu.

R0 Basic reproductive rate -
The expected number of
secondary cases produced by
a single (typical) infection
in a completely susceptible
population. It reflects the
ability to transmit the
disease.
R0 tells you the average
number of people who will
contract a contagious
disease from one person with
that disease. It
specifically applies to a
population of people who
were previously free of
infection and haven’t been
vaccinated.

R0 Seasonal influenza = 1-2,
COVID-19 = 2-3, SARS = 3-5,
Measles = 12-18

if value of R0:
<1: the transmission will
wane off because one
infectious case will infect
less than one person on
average.
>1: the transmission is
likely to continue in a
population

R0 is based on three
factors: 
R0 = T * C * D

- Transmissibility (T) is
the probability of infection
given contact between a
susceptible and infected
individual Any easily
transmissible infection will
have a higher R0. For
example, an airborne
infection like measles, flu
has a higher R0 than others
requiring direct contacts
like HIV or Ebola.
- Average rate of contact
(C) between susceptible and
infected individuals
- Duration (D) of
infectiousness.

Virulenceisa pathogen's
ability to cause disease. It
determines the degree of
pathogenicity.
The infectious period refers
to the duration of time
during which an individual
infected with a contagious
disease, such as the flu,
can transmit the infection
to others. A long period of
infectiousness will
contribute to a higher R0
value.
The contact rate, also known
as the contact frequency or
contact rate parameter,
refers to the average number
of contacts or interactions
that an individual has with
others within a specific
period, typically within a
day. A high contact rate
will contribute to a higher
R0 value.
Mode of transmission: method
of transmission (direct
contact, airborne, droplet,
fecal-oral))");

const ScriptTemplate * ScriptTemplate::Epidemiology() {return &EpidemiologyScriptTemplate;}

constexpr ScriptTemplate ExportSynthesizedProteinScriptTemplate("ExportSynthesizedProtein.py", "\x01" R"(Once the protein is
synthesized, it is
transported to the
ENDOPLASMIC RETICULUM (ER),
where it undergoes further
modifications, such as the
addition of carbohydrate
side chains. The ER also
makes phospholipids for
other cellular membranes,
which are transported when
the VESICLE forms. If the
modified proteins are not
destined to stay in the ER,
they will be packaged into
vesicles, or small spheres
of membrane that are used
for transport, and shipped
to the GOLGI APPARATUS. The
GOLGI APPARATUS functions as
a molecular assembly line in
which membrane proteins
undergo extensive
post-translational
modification.
Finally, the protein is
transported to its final
destination, which could be
the cell membrane,
lysosomes, or outside the
cell. The CYTOSKELETON
provides a network of
protein fibers that help to
transfer VESICLES to their
proper locations.

)");

const ScriptTemplate * ScriptTemplate::ExportSynthesizedProtein() {return &ExportSynthesizedProteinScriptTemplate;}

constexpr ScriptTemplate FluScriptTemplate("Flu.py", "\x01" R"(Influenza, commonly known as
the flu, is caused by
influenza viruses. These
viruses belong to the
Orthomyxoviridae family and
are categorized into types
A, B, C, and D. Types A and
B are responsible for
seasonal flu outbreaks in
humans, while type C usually
causes mild respiratory
illness and type D primarily
affects cattle.
The influenza virus is
unique because of its
ability to rapidly evolve
through mutations and
genetic reassortment,
leading to the emergence of
new strains. This process,
known as antigenic drift and
antigenic shift,
respectively, is why
seasonal flu vaccines need
to be updated regularly to
match the circulating
strains.
The flu primarily spreads
between humans through
respiratory droplets
produced when an infected
person coughs, sneezes, or
talks. These droplets can
land in the mouths or noses
of people who are nearby
(usually within about 6
feet), or they can be
inhaled into the lungs. This
is considered the main mode
of transmission for
influenza viruses.
Additionally, people can
also become infected with
the flu virus by touching
surfaces or objects that
have the virus on them and
then touching their own
mouth, nose, or possibly
their eyes.
Influenza viruses primarily
target the respiratory
tract, infecting epithelial
cells in the nose, throat,
and lungs. The virus spreads
through respiratory droplets
when an infected person
coughs, sneezes, or talks.
It can also survive on
surfaces for a short period,
allowing indirect
transmission through contact
with contaminated objects.
Symptoms of influenza
typically include fever,
cough, sore throat, muscle
aches, fatigue, and
sometimes gastrointestinal
symptoms like nausea and
vomiting. While most people
recover from the flu without
complications, it can lead
to severe illness,
particularly in young
children, older adults,
pregnant women, and
individuals with underlying
health conditions.
Preventive measures such as
vaccination, good hygiene
practices (like frequent
handwashing), and avoiding
close contact with sick
individuals can help reduce
the spread of the flu virus
between humans and animals.
Additionally, antiviral
medications may be
prescribed to treat the flu
and lessen its severity if
taken early in the illness.
Good Hygiene Practices, like
frequent hand washing, using
of sanitizers, avoid
touching face, eyes, nose,
mouth, and cover the mouth
and nose when coughing or
sneezing, help to reduce the
spread of the flu.
Treating the flu involves a
multifaceted approach to
alleviate symptoms and
support the body's immune
response. Rest and hydration
are paramount, allowing the
body to conserve energy and
combat the infection
effectively.
Over-the-counter medications
like pain relievers and
decongestants can provide
relief from fever, body
aches, and nasal congestion,
while antiviral medications
may be prescribed in severe
cases to reduce the duration
and severity of the illness.
Supportive measures such as
steam inhalation, nasal
saline irrigation, and
herbal remedies can also
help alleviate symptoms and
promote comfort during
recovery.)");

const ScriptTemplate * ScriptTemplate::Flu() {return &FluScriptTemplate;}

constexpr ScriptTemplate GeneMutationsScriptTemplate("GeneMutations.py", "\x01" R"(Cells containing an abnormal
chromosome structure within
a gene due to a mistake in
one (point mutation) or more
pairs of nucleotides. It can
lead to the production of an
abnormal protein.)");

const ScriptTemplate * ScriptTemplate::GeneMutations() {return &GeneMutationsScriptTemplate;}

constexpr ScriptTemplate GeneRegulationScriptTemplate("GeneRegulation.py", "\x01" R"(Prokaryotes and eukaryotes
alter gene expression in
response to their changing
environment. In
multicellular eukaryotes,
gene expression regulates
development and is
responsible for differences
in cell types. RNA molecules
play many roles in
regulating gene expression
in eukaryotes.
Natural selection has
favored bacteria that
produce only the products
needed by that cell. A cell
can regulate the production
of enzymes by feedback
inhibition or by gene
regulation. Gene expression
in bacteria is controlled by
the operon model.)");

const ScriptTemplate * ScriptTemplate::GeneRegulation() {return &GeneRegulationScriptTemplate;}

constexpr ScriptTemplate GenomeMutationsScriptTemplate("GenomeMutations.py", "\x01" R"(Inpolyploidy all
pairsofchromosomes are
affected and in aneuploidy
only affects 1 pair of
chromosomes is affected.
A KARYOGRAM is a graphical
representation of
a KARYOTYPE, which is the
complete set of chromosomes
in an organism or cell. The
chromosomes are generally
organized in pairs, ordered
by size and position of
centromere for chromosomes
of the same size.
Karyotyping generally
combines light microscopy
and photography in the
metaphase of the cell cycle,
and results in a
photomicrographic karyogram.
In contrast, a schematic
KARYOGRAM is a designed
graphic representation of a
KARYOTYPE.
KARYOTYPE FOLRMULAS:
Monosomy:
45, X (Turner Syndrome)
Trisomy:
47, XXY (Klinefelter
syndrome)
47, XX +21 (Down syndrome
female)
47, XY +21 (Down syndrome
male))");

const ScriptTemplate * ScriptTemplate::GenomeMutations() {return &GenomeMutationsScriptTemplate;}

constexpr ScriptTemplate HelperTCellsScriptTemplate("HelperTCells.py", "\x01" R"(Helper T cells are a
critical component of the
adaptive immune system
responsible for
orchestrating and regulating
immune responses. They play
a central role in
coordinating the activities
of other immune cells,
including B cells, cytotoxic
T cells, and macrophages, to
mount effective immune
responses against pathogens.
Functions:
1. Activation of Immune
Responses: Helper T cells
are activated when they
encounter antigens presented
by antigen-presenting cells
(APCs) such as dendritic
cells, macrophages, or B
cells. This activation
usually occurs in lymphoid
organs like lymph nodes or
the spleen. Upon activation,
helper T cells undergo
clonal expansion,
proliferating to generate a
large population of effector
cells.
2. Cytokine Production:
Activated helper T cells
produce a variety of
cytokines, signaling
molecules that regulate the
immune response. Different
subsets of helper T cells
produce distinct cytokine
profiles that influence the
nature of the immune
response.
3. Activation of B Cells:
Helper T cells provide
crucial help to B cells
during the antibody
response. They stimulate B
cells to undergo
proliferation and
differentiation into plasma
cells, which secrete
antibodies specific to the
antigen. Additionally,
helper T cells promote class
switching in B cells,
resulting in the production
of different antibody
isotypes with distinct
effector functions.
4. Activation of Cytotoxic T
Cells: Helper T cells
facilitate the activation
and differentiation of
cytotoxic T cells, which are
responsible for directly
killing infected cells.
Helper T cells provide
signals that promote the
expansion and
differentiation of cytotoxic
T cells into effector cells
capable of recognizing and
eliminating cells harboring
intracellular pathogens.
5. Regulation of Immune
Responses: Helper T cells
also play a role in
regulating immune responses
to prevent excessive
inflammation and tissue
damage. Regulatory T cells
suppress immune responses
and maintain immune
tolerance by inhibiting the
activation and function of
other immune cells.)");

const ScriptTemplate * ScriptTemplate::HelperTCells() {return &HelperTCellsScriptTemplate;}

constexpr ScriptTemplate HeritabilityScriptTemplate("Heritability.py", "\x01" R"(Heritability is a concept in
genetics that quantifies the
extent to which genetic
variation contributes to
phenotypic variation within
a population. It represents
the proportion of phenotypic
variation that can be
attributed to genetic
differences among
individuals, expressed as a
percentage.
To deduce the concept of
heritability, consider a
population where a
particular trait, such as
height, is known to vary
among individuals.
Heritability would describe
the extent to which
differences in genetic
makeup contribute to the
observed variation in height
within that population.
For example, if the
heritability of height is
80%, it means that 80% of
the observed variation in
height among individuals in
the population is due to
genetic differences, while
the remaining 20% is
influenced by environmental
factors, random
developmental processes, and
measurement error.)");

const ScriptTemplate * ScriptTemplate::Heritability() {return &HeritabilityScriptTemplate;}

constexpr ScriptTemplate HIVScriptTemplate("HIV.py", "\x01" R"(Human immunodeficiency virus
infects helper T cells. The
loss of helper T cells
impairs both the humoral and
cell mediated immune
responses and leads to AIDS.
HIV eludes the immune system
because of antigenic
variation and an ability to
remain latent while
integrated into host DNA.
People with AIDS are highly
susceptible to opportunistic
infections and cancers that
take advantage of an immune
system in collapse. The
spread of HIV is a worldwide
problem The best approach
for slowing this spread is
education about practices
that transmit the virus.)");

const ScriptTemplate * ScriptTemplate::HIV() {return &HIVScriptTemplate;}

constexpr ScriptTemplate HumoralResponseScriptTemplate("HumoralResponse.py", "\x01" R"(The humoral response is a
critical component of the
body's immune system
responsible for defending
against pathogens,
particularly those that are
circulating freely in bodily
fluids like blood and lymph.
This aspect of immunity
involves the production of
specific proteins called
antibodies, which are
produced by a type of white
blood cell called B
lymphocytes or B cells.
During an immune response,
when the body encounters a
foreign antigen (a molecule
recognized as foreign by the
immune system), B cells are
activated and undergo a
process called clonal
expansion. This process
involves the proliferation
and differentiation of B
cells into plasma cells,
which are antibody-secreting
cells.
Plasma cells produce large
quantities of antibodies
that are specifically
targeted against the antigen
encountered. These
antibodies circulate in the
bloodstream and other bodily
fluids, binding to and
neutralizing pathogens such
as bacteria, viruses, and
toxins. Antibodies can
prevent pathogens from
infecting host cells, mark
them for destruction by
other immune cells, or
facilitate their removal by
other parts of the immune
system.
The humoral response is
particularly effective
against extracellular
pathogens that are present
in bodily fluids or tissues
outside of cells. It plays a
crucial role in protecting
against bacterial
infections, some viral
infections, and certain
toxins. Additionally,
humoral immunity provides
long-term protection against
pathogens through the
production of memory B
cells, which can quickly
mount a secondary immune
response upon re-exposure to
the same antigen, leading to
a more rapid and robust
antibody production.
Vaccines work by stimulating
the humoral immune response,
inducing the production of
specific antibodies against
the targeted pathogen.
Understanding the humoral
response is essential for
the development of vaccines,
immunotherapies, and
treatments for infectious
diseases and autoimmune
disorders.)");

const ScriptTemplate * ScriptTemplate::HumoralResponse() {return &HumoralResponseScriptTemplate;}

constexpr ScriptTemplate ImmuneRejectionScriptTemplate("ImmuneRejection.py", "\x01" R"(Cells transferred from one
person to another can be
attacked by immune defenses.
This complicates blood
transfusions or the
transplant of tissues or
organs.

Blood group:
Antigens on red blood cells
determine whether a person
has blood type A (A
antigen), B (B antigen), AB
(both A and B antigens), or
O (neither antigen).
Antibodies to nonself blood
types exist in the body.
Transfusion with
incompatible blood leads to
destruction of the
transfused cells. Recipient
donor combinations can be
fatal or safe.

Tissue and organ
transplants: 
MHC molecules are different
among genetically
nonidentical individuals.
Differences in MHC molecules
stimulate rejection of
tissue grafts and organ
transplants. Chances of
successful transplantation
increase if donor and
recipient MHC tissue types
are well matched.
Immunosuppressive drugs
facilitate transplantation.
Lymphocytes in bone marrow
transplants may cause the
donor tissue to reject the
recipient.)");

const ScriptTemplate * ScriptTemplate::ImmuneRejection() {return &ImmuneRejectionScriptTemplate;}

constexpr ScriptTemplate ImmunologicalMemoryScriptTemplate("ImmunologicalMemory.py", "\x01" R"(Immunological memory is
responsible for long term
protections against
diseases, due to either a
prior infection or
vaccination. The first
exposure to a specific
antigen represents the
primary immune response.
During this time, selected B
and T cells give rise to
their effector forms. In the
secondary immune response,
memory cells facilitate a
faster, more efficient
response.)");

const ScriptTemplate * ScriptTemplate::ImmunologicalMemory() {return &ImmunologicalMemoryScriptTemplate;}

constexpr ScriptTemplate InflammatoryResponsesScriptTemplate("InflammatoryResponses.py", "\x01" R"(The inflammatory response ,
such as pain and swelling,
is brought about by
molecules released upon
injury of infection.
Mast cells, a type of
connective tissue, release
histamine , which triggers
blood vessels to dilate and
become more permeable.
Activated macrophages and
neutrophils release
cytokines , signaling
molecules that enhance the
immune response.
Pus, a fluid rich in white
blood cells, dead pathogens,
and cell debris from damaged
tissues.
Inflammation can be either
local or systemic
(throughout the body)
examples of inflammatory
responses:
- Fever is a systemic
inflammatory response
triggered by pyrogens
released by macrophages and
by toxins from pathogens. 
- Septic shock is a life
threatening condition caused
by an overwhelming
inflammatory response.
- Swelling: isolates dead
cells, attraction of immune
cells, prevention of future
damage, swollen tissue
begins to shrink to reduce
blood flow to the area.
- Redness: more blood to the
affected area (more oxygen,
nutrients, and immune cells)
- Heat: Enhanced metabolic
activity, acceleration of
healing
- Pain: warning sign,
limiting movement,
initiation of protective
reflexes)");

const ScriptTemplate * ScriptTemplate::InflammatoryResponses() {return &InflammatoryResponsesScriptTemplate;}

constexpr ScriptTemplate InnateImmunityOfInvertebratesScriptTemplate("InnateImmunityOfInvertebrates.py", "\x01" R"(In insects, an exoskeleton
made of chitin forms the
first barrier to pathogens.
The digestive system is
protected by a chitin based
barrier and lysozyme, an
enzyme that breaks down
bacterial cell walls.
Hemocytes (a kind of a cell)
circulate within hemolymph
and carry out phagocytosis,
the ingestion and digestion
of foreign substances
including bacteria.
Hemocytes also secrete
antimicrobial peptides that
disrupt the plasma membranes
of fungi and bacteria.)");

const ScriptTemplate * ScriptTemplate::InnateImmunityOfInvertebrates() {return &InnateImmunityOfInvertebratesScriptTemplate;}

constexpr ScriptTemplate InnateImmunityOfVertebratesScriptTemplate("InnateImmunityOfVertebrates.py", "\x01" R"(see Barrier...,
CellularInnate...,
Antimicrobal...,
Inflammatory..., Adaptive...

Some pathogens avoid
destruction by modifying
their surface to prevent
recognition or by resisting
breakdown following
phagocytosis. Tuberculosis
is one such disease and
kills more than a million
people a year)");

const ScriptTemplate * ScriptTemplate::InnateImmunityOfVertebrates() {return &InnateImmunityOfVertebratesScriptTemplate;}

constexpr ScriptTemplate InterleukinsScriptTemplate("Interleukins.py", "\x01" R"(Helper T cells, are capable
of producing interleukins,
which are signaling
molecules crucial for
orchestrating various
aspects of the adaptive
immune response.
Interleukins play diverse
roles in modulating the
functions of immune cells,
coordinating their
activities, and shaping the
overall immune response.
Interleukins are signaling
molecules that promote a
global adapted immune
response.)");

const ScriptTemplate * ScriptTemplate::Interleukins() {return &InterleukinsScriptTemplate;}

constexpr ScriptTemplate LactoseOperonScriptTemplate("LactoseOperon.py", "\x01" R"(The operon regulates and
expresses
enzymes/transporters for the
synthesis/uptake/metabolise
the disaccharide lactose
into galactose in
prokaryotes.
Repressor status when
synthesised by the regulator
gene: Lactose is not usually
available on the
extracellular region.
Prokaryotes don’t need to
metabolise/uptake it (they
rather use glucose!). The
repressor is usually
synthesised ACTIVE.
Corepressor/Inducer (What is
the name? What type of
substance?):
INDUCER Allolactose
Is the operon usually ON or
OFF?: OFF
Status change What condition
triggers a status change?:
When there is no glucose and
only lactose!
Corepressor/Inducer
role Where does the
corepressor/Inducer binds?
What is causes?: Allolactose
binds to the allosteric site
of the repressor changing
its conformation into the
INACTIVE form.)");

const ScriptTemplate * ScriptTemplate::LactoseOperon() {return &LactoseOperonScriptTemplate;}

constexpr ScriptTemplate LymphaticSystemScriptTemplate("LymphaticSystem.py", "\x01" R"(Memory B-Cells:
- After encountering
antigens and undergoing
activation in secondary
lymphoid organs such as
lymph nodes or the spleen, a
subset of B-cells
differentiates into memory
B-cells.
- Memory B-cells can then
migrate to and reside within
various lymphoid tissues,
including lymph nodes, the
spleen, and
mucosal-associated lymphoid
tissues (MALT), such as
tonsils and Peyer's patches
in the intestines.
- Within these lymphoid
tissues, memory B-cells
persist for extended
periods, ready to rapidly
respond to re-exposure to
the same antigen by quickly
differentiating into
antibody-secreting plasma
cells, thereby providing an
accelerated and more robust
immune response during
subsequent encounters with
the pathogen.
Memory T-Cells:
- Similarly, memory T-cells
are generated following
activation and
differentiation in secondary
lymphoid organs,
particularly the thymus for
T-cells.
- Memory T-cells can migrate
and distribute throughout
the body via the bloodstream
and lymphatic vessels,
allowing them to survey
various tissues and organs
for the presence of
pathogens or infected cells.
- Like memory B-cells,
memory T-cells can reside
within lymphoid tissues, as
well as peripheral tissues
where they may encounter
pathogens.
- Upon re-exposure to the
same antigen, memory T-cells
can rapidly respond by
proliferating and
differentiating into
effector T-cells, which
exert their immune
functions, such as cytokine
production or direct killing
of infected cells, thus
providing enhanced
protection against the
pathogen.)");

const ScriptTemplate * ScriptTemplate::LymphaticSystem() {return &LymphaticSystemScriptTemplate;}

constexpr ScriptTemplate MechanismusOfPostTranscriptionalRegulationScriptTemplate("MechanismusOfPostTranscriptionalRegulation.py", "\x01" R"(Transcription alone does not
account for gene expression.
Regulatory mechanisms can
operate at various stages
after transcription. Such
mechanisms allow a cell to
fine-tune gene expression
rapidly in response to
environmental changes.

RNA Processing: In
alternative RNA splicing,
different mRNA molecules are
produced from the same
primary transcript,
depending on which RNA
segments are treated as
exons and which as introns.

mRNA Degradation: The life
span of mRNA molecules in
the cytoplasm is a key to
determining protein
synthesis. Eukaryotic mRNA
is more long lived than
prokaryotic mRNA. Nucleotide
sequences that influence the
lifespan of mRNA in
eukaryotes reside in the
untranslated region (UTR) at
the 3end of the molecule.

Initiation of Translation:
The initiation of
translation of selected
mRNAs can be blocked by
regulatory proteins that
bind to sequences or
structures of the mRNA.
Alternatively, translation
of all mRNAs in a cell may
be regulated simultaneously.
For example, translation
initiation factors are
simultaneously activated in
an egg following
fertilization.

Protein Processing and
Degradation: After
translation, various types
of protein processing,
including cleavage and the
addition of chemical groups,
are subject to control.
Proteasomes are giant
protein complexes that bind
protein molecules and
degrade them.)");

const ScriptTemplate * ScriptTemplate::MechanismusOfPostTranscriptionalRegulation() {return &MechanismusOfPostTranscriptionalRegulationScriptTemplate;}

constexpr ScriptTemplate MemoryBCellsScriptTemplate("MemoryBCells.py", "\x01" R"(Memory B cells are a crucial
component of the adaptive
immune system, playing a key
role in the establishment of
immunological memory and
long-term protection against
pathogens. These cells are a
subset of B lymphocytes, a
type of white blood cell,
and are generated during an
immune response to an
antigen.
When the immune system
encounters a foreign
antigen, such as a virus or
bacteria, B cells are
activated and undergo
proliferation and
differentiation. Some of
these activated B cells
differentiate into plasma
cells, which produce and
secrete antibodies to
neutralize the antigen.
However, a portion of the
activated B cells
differentiate into memory B
cells instead of plasma
cells.
Сharacteristics:
1. Longevity: Memory B cells
can persist in the body for
months to years after the
initial infection or
immunization, providing
long-term immunity against
specific pathogens.
2. Quiescence: Unlike plasma
cells, which are actively
secreting antibodies, memory
B cells are quiescent and do
not produce antibodies
unless they are reactivated
by exposure to the same
antigen.
3. Rapid Response: Upon
re-exposure to the antigen,
memory B cells can quickly
become activated and
differentiate into
antibody-secreting plasma
cells. This rapid response
enables the immune system to
mount a faster and more
robust defense against
recurrent infections.
4. Affinity Maturation:
Memory B cells often undergo
affinity maturation during
the initial immune response,
resulting in the production
of antibodies with higher
affinity for the antigen.
This process contributes to
the effectiveness of the
secondary immune response
upon re-exposure to the
antigen.
Immunization stimulates the
generation of memory B cells
specific to the antigens
present in the vaccine. In
the event of subsequent
exposure to the pathogen,
these memory B cells can
rapidly recognize and
respond to the antigen,
leading to a quicker and
more effective immune
response that can prevent or
mitigate the severity of the
infection.)");

const ScriptTemplate * ScriptTemplate::MemoryBCells() {return &MemoryBCellsScriptTemplate;}

constexpr ScriptTemplate MemoryTCellsScriptTemplate("MemoryTCells.py", "\x01" R"(Memory T-Cells play a
crucial role in adaptive
immunity by providing
long-lasting protection
against recurrent infections
and contributing to the
rapid and robust response
upon re-exposure to the same
antigen.
Functions:
1. Formation: Memory T cells
are generated during the
primary immune response
following exposure to an
antigen. After activation by
antigen-presenting cells
(APCs), such as dendritic
cells, naive T cells
differentiate into effector
T cells, which carry out the
immediate response against
the antigen, and a subset of
these effector T cells
further differentiate into
memory T cells.
2. Longevity: Memory T cells
can persist in the body for
extended periods, ranging
from months to years,
following the resolution of
the initial infection or
immunization. This long-term
persistence allows memory T
cells to provide durable
protection against
reinfection with the same
pathogen.
3. Quiescence and Readiness:
Unlike effector T cells,
which are actively involved
in the immune response,
memory T cells are quiescent
and do not exhibit effector
functions under normal
conditions. However, they
remain poised and "ready to
go," quickly responding to
antigen re-encounter by
rapidly proliferating and
differentiating into
effector cells.
4. Rapid Response: Memory T
cells mount a faster and
more robust immune response
upon re-exposure to the
antigen compared to naive T
cells. This rapid response
is facilitated by several
factors, including the
presence of pre-formed
effector molecules, enhanced
sensitivity to antigen
stimulation, and increased
expression of co-stimulatory
molecules.
5. Heterogeneity: Memory T
cells can be broadly
categorized into two main
subsets based on the
expression of specific
surface markers and
functional properties:
central memory T cells (Tcm)
and effector memory T cells
(Tem). Central memory T
cells primarily reside in
secondary lymphoid organs
and exhibit high
proliferative capacity upon
antigen re-encounter,
whereas effector memory T
cells are found in
peripheral tissues and
provide immediate effector
functions at sites of
infection.
6. Immunological Memory: The
presence of memory T cells
enables the immune system to
establish immunological
memory, which is the ability
to mount a heightened and
more effective immune
response upon subsequent
exposures to the same
antigen. This phenomenon
underlies the basis of
vaccination, where the
introduction of antigenic
material primes the immune
system to generate memory T
cells and confer long-term
protection against specific
pathogens.)");

const ScriptTemplate * ScriptTemplate::MemoryTCells() {return &MemoryTCellsScriptTemplate;}

constexpr ScriptTemplate MHCIScriptTemplate("MHCI.py", "\x01" R"(The major histocompatibility
complex (MHC) class I cell
surface protein complex
plays a fundamental role in
the immune system by
presenting peptide fragments
derived from intracellular
proteins to cytotoxic T
lymphocytes (CTLs). This
process is crucial for the
recognition of self by
leukocytes and the
maintenance of immune
surveillance against
infected or abnormal cells.
The role of MHC class I in
recognition of self by
leukocytes:
1. Intracellular Protein
Processing: Within cells,
proteins are continuously
degraded into peptide
fragments by proteasomes,
specialized
protein-degrading complexes.
These peptide fragments are
typically derived from
endogenous proteins
synthesized within the cell,
including self-proteins and
foreign proteins, such as
those from viruses or
intracellular bacteria.
2. Peptide Loading onto MHC
Class I Molecules: Peptide
fragments generated by
proteasomal degradation are
transported into the
endoplasmic reticulum (ER)
by transporter proteins.
Within the ER, these peptide
fragments are loaded onto
newly synthesized MHC class
I molecules, which consist
of a peptide-binding groove
formed by the MHC heavy
chain and a beta-2
microglobulin protein.
3. Surface Expression of MHC
Class I-Peptide Complexes:
Peptide-loaded MHC class I
molecules are transported
from the ER to the cell
surface via the Golgi
apparatus. Once at the cell
surface, these MHC class
I-peptide complexes are
displayed on the plasma
membrane, where they can be
recognized by immune cells,
particularly cytotoxic T
lymphocytes (CTLs).
4. Recognition by Cytotoxic
T Lymphocytes (CTLs): CTLs
express T cell receptors
(TCRs) that can specifically
recognize peptide antigens
presented by MHC class I
molecules on target cells.
When a CTL encounters a cell
displaying foreign peptide
fragments (e.g., from a
virus-infected cell or a
tumor cell), its TCR binds
to the MHC class I-peptide
complex, triggering the CTL
to initiate its effector
functions, such as releasing
cytotoxic molecules that
induce apoptosis in the
target cell.
5. Immune Surveillance and
Self-Recognition: By
presenting peptide fragments
derived from endogenous
proteins, including
self-proteins, on their
surface, cells express a
"snapshot" of their internal
protein content to CTLs.
Healthy cells typically
present self-peptides
derived from normal cellular
proteins, allowing CTLs to
distinguish between self and
non-self. Abnormal or
infected cells, however, may
present foreign or aberrant
peptides derived from viral
or mutated proteins,
triggering CTL-mediated
immune responses against
these cells.)");

const ScriptTemplate * ScriptTemplate::MHCI() {return &MHCIScriptTemplate;}

constexpr ScriptTemplate MHCIIScriptTemplate("MHCII.py", "\x01" R"(1. Expression: MHC II
molecules are constitutively
expressed on the surface of
professional
antigen-presenting cells
(APCs), such as dendritic
cells, macrophages, and B
cells. Expression can be
upregulated by various
stimuli, including
inflammatory signals and
cytokines.
2. Antigen Presentation: MHC
II molecules primarily
present peptide fragments
derived from exogenous
antigens that have been
internalized by the APC
through processes such as
phagocytosis, endocytosis,
or pinocytosis. Within the
endosomal compartments of
the APC, these antigens are
proteolytically degraded
into smaller peptides, which
then bind to the
peptide-binding groove of
MHC II molecules.
3. Interaction with Helper T
Cells: Once loaded with
antigenic peptides, MHC II
molecules are transported to
the cell surface, where they
present these peptides to
CD4+ helper T cells. The T
cell receptor (TCR) on the
surface of the CD4+ T cell
recognizes the peptide-MHC
II complex, leading to T
cell activation.
4. Helper T Cell Activation:
The interaction between the
TCR on the CD4+ T cell and
the peptide-MHC II complex
on the APC, along with
co-stimulatory signals
provided by molecules on the
APC, leads to the activation
of the helper T cell.
Activated helper T cells
then secrete cytokines that
regulate and orchestrate
immune responses, such as
promoting B cell activation,
enhancing macrophage
function, and recruiting
other immune cells to the
site of infection.
5. Role in Immune Responses:
MHC II-mediated antigen
presentation is essential
for the initiation and
regulation of adaptive
immune responses against
extracellular pathogens,
including bacteria,
parasites, and fungi. It
enables the immune system to
recognize and respond to a
wide range of foreign
antigens, contributing to
host defense and immune
surveillance.)");

const ScriptTemplate * ScriptTemplate::MHCII() {return &MHCIIScriptTemplate;}

constexpr ScriptTemplate NoncodingRNAInControllingGeneExpressionScriptTemplate("NoncodingRNAInControllingGeneExpression.py", "\x01" R"(Only a small fraction of DNA
codes for proteins, and a
very small fraction of the
non-protein-coding DNA
consists of genes for RNA
such as rRNA and tRNA. A
significant amount of the
genome may be transcribed
into noncoding RNAs
(ncRNAs). Noncoding RNAs
regulate gene expression at
two points: mRNA translation
and chromatin configuration.

MicroRNAs (miRNAs) are small
single-stranded RNA
molecules that can bind to
mRNA. These can degrade mRNA
or block its translation.
The phenomenon of inhibition
of gene expression by RNA
molecules is called RNA
interference(RNAi). RNAi is
caused by small interfering
RNAs (siRNAs). siRNAs and
miRNAs are similar but form
from different RNA
precursors.

Small ncRNAs can regulate
gene expression at multiple
steps. An increase in the
number of miRNAs in a
species may have allowed
morphological complexity to
increase over evolutionary
time. siRNAs may have
evolved first, followed by
miRNAs and later piRNAs.)");

const ScriptTemplate * ScriptTemplate::NoncodingRNAInControllingGeneExpression() {return &NoncodingRNAInControllingGeneExpressionScriptTemplate;}

constexpr ScriptTemplate OperonsScriptTemplate("Operons.py", "\x01" R"(A cluster of functionally
related genes can be under
coordinated control by a
single “on-off switch”. The
regulatory “switch” is a
segment of DNA called an
operator usually positioned
within the promoter. An
operon is the entire stretch
of DNA that includes the
operator, the promoter, and
the genes that they control.
The operon can be switched
off by a protein repressor.
The repressor prevents gene
transcription by binding to
the operator and blocking
RNA polymerase. The
repressor is the product of
a separate regulatory gene.
The repressor can be in an
active or inactive form,
depending on the presence of
other molecules. A
corepressoris a molecule
that cooperates with a
repressor protein to switch
an operon off. For example,
E. coli can synthesize the
amino acid tryptophan.
A repressible operon is one
that is usually on; binding
of a repressor to the
operator shuts off
transcription. The trp
operon is a repressible
operon. An inducible operon
is one that is usually off;
a molecule called an inducer
inactivates the repressor
and turns on transcription.
The lacoperon is an
inducible operon and
contains genes that code for
enzymes used in the
hydrolysis and metabolism of
lactose. By itself, the lac
repressor is active and
switches the lac operon off.
A molecule called an inducer
inactivates the repressor to
turn the lac operon on.
Inducible enzymes usually
function in catabolic
pathways; their synthesis is
induced by a chemical
signal. Repressible enzymes
usually function in anabolic
pathways; their synthesis is
repressed by high levels of
the end product. Regulation
of the trpand lacoperons
involves negative control of
genes because operons are
switched off by the active
form of the repressor.)");

const ScriptTemplate * ScriptTemplate::Operons() {return &OperonsScriptTemplate;}

constexpr ScriptTemplate OxidizedReducedScriptTemplate("OxidizedReduced.py", "\x01" R"(Reduced molecule loses
electrons
Oxidized molecule gains
electrons)");

const ScriptTemplate * ScriptTemplate::OxidizedReduced() {return &OxidizedReducedScriptTemplate;}

constexpr ScriptTemplate PenetranceScriptTemplate("Penetrance.py", "\x01" R"(Penetrance in genetics
refers to the proportion of
individuals with a specific
genotype (genetic makeup)
who exhibit the associated
phenotype (observable traits
or characteristics) in a
population. It is a
statistical measure that
quantifies the likelihood or
probability of a particular
genotype manifesting as the
corresponding phenotype.)");

const ScriptTemplate * ScriptTemplate::Penetrance() {return &PenetranceScriptTemplate;}

constexpr ScriptTemplate PhagocytosisScriptTemplate("Phagocytosis.py", "\x01" R"(Phagocytosis is a process by
which a cell uses its plasma
membrane to engulf a
particle, giving rise to an
internal compartment called
the phagosome. It is one
type of endocytosis. A cell
that performs phagocytosis
is called a phagocyte.
Cells of the immune system
can perform phagocytosis:
neutrophils, macrophages,
dendritic cells, and B
lymphocytes. The act of
phagocytizing pathogenic or
foreign particles allows
cells of the immune system
to know what they are
fighting against.
Steps:
1. Activation and
chemotaxis: Phagocytes are
attracted to the site of
infection by various
substances generated in the
immune response. Resting
phagocytes are activated by
inflammatory mediators,
which increase their
metabolic and microbicidal
activity. Activated cells
also express more
glycoprotein receptors,
which help them to reach the
site of infections as well
as to bind firmly with
microorganisms.
2. Recognition of invading
microbes: The next step in
phagocytosis is the
adherence of the antigen to
the cell membrane of the
phagocytic cells. Adherence
induces membrane
protrusions, called
pseudopodia, to extend
around the attached material
and to ingest them.
Phagocytic cells contain
various receptors which help
them to attach with bacteria
or viruses.
3. Ingestion and formation
of phagosomes: The
phagocytic cell engulfs the
microbe or particle, forming
a phagosome. The phagosome
is a vesicle that contains
the ingested material.
4. Formation of phagolysome:
The phagosome fuses with
lysosomes, forming a
phagolysosome. The lysosomes
contain digestive enzymes
that break down the ingested
material.
5. Microbial killing and
formation of residual
bodies: The ingested
material is digested by the
lysosomal enzymes, and the
residual bodies are formed.
The residual bodies are then
exocytosed from the cell
with bound MHC II.)");

const ScriptTemplate * ScriptTemplate::Phagocytosis() {return &PhagocytosisScriptTemplate;}

constexpr ScriptTemplate PhotosynthesisScriptTemplate("Photosynthesis.py", "\x01" R"(Formula: 6CO2 + 6H2O ->
(with sunlight and
chlorophyll) C6H12O6 + 6O2
The amount of ATP produced
during photosynthesis is not
a fixed value, as it depends
on various factors such as
the intensity of light, the
type of plant, and the
availability of water and
carbon dioxide.
The chemical process by
which plants use
CARBON-DIOXIDE (CO2), WATER
and ENERGY (SUNLIGHT) to
manufacture GLUCOSE, the
building blocks of plants,
is called photosynthesis. In
process, OXYGEN is produced
as a byproduct.
Photosynthesis occurs in
CHLOROPLASTS.
Photosynthesis consists of
two sets of reactions - the
LIGHT DEPENDENT REACTIONS
and the CALVIN CYCLE
(LIGHT-INDEPENDENT
REACTIONS).
The Calvin cycle occurs in
the STROMA, the
light-dependent reactions
occur in the THYLAKOID.
Thylakoids contains pairs of
PHOTOSYSTEMS, called
photosystem 1 and
photosystem 2, that work in
tandem to produce the energy
that will later be used in
the stroma to manufacture
sugars. 
The photosystems of the
thylakoid consist of a
network of accessory pigment
molecules and chlorophyll -
the molecule that absorb the
photons of light. The gained
light energy excites
electrons to a higher state.
Photosystems will direct the
excitation energy gathered
by the pigment molecules to
a reaction center chlorophyl
molecule which will then
pass the electrons to a
series of proteins located
on the thylakoid membrane.
#### In PS 2
The energized electrons are
passed from the reaction
center of a PS 2 to an
electron transport chain.
The electrons lost by PS 2
are replaced by a process
called photolysis, which
involves the oxidation of a
water molecule, producing
free electrons and oxygen
gas. 
As electrons pass through
the electron transport
chain, the energy from the
electron is used to pump
hydrogen ions (H+) from the
stroma to the thylakoid,
creating a CONCENTRATION
GRADIENT. This gradient
powers a protein called
ATP-SYNTHASE, which
phosphorylates ADP to form
ATP.
#### In PS 1
The low-energy electrons
leaving PS 2 are shuttled to
PS 1. There, low energy
electrons are reenergized
and are passed through an
electron transport chain
where they are used to
reduce the electron carrier
NADP+ to NADPH.
NADPH and ATP are provided
to the metabolic pathways in
the stroma. Therefore these
ATP and NADPH are used in
the stroma to fuel the
Calvin Cycle reactions.
The Calvin Cycle consists of
a series of reactions that
reduce carbon dioxide to
produce the carbohydrate
glyceraldehyde-3-phosphate
(G3P). 
The Calvin cycle has to run
6 times to produce one
molecule of glucose. These
molecules can remove their
phosphate and add fructose
to form sucrose, the
molecule plants use to
transport carbohydrates
throughout their system.
Glucose phosphate is also
the starting molecule for
the synthesis of starch and
cellulose.

)");

const ScriptTemplate * ScriptTemplate::Photosynthesis() {return &PhotosynthesisScriptTemplate;}

constexpr ScriptTemplate PlasmaCellsScriptTemplate("PlasmaCells.py", "\x01" R"(Plasma cells are a
specialized type of white
blood cell, or lymphocyte,
that plays a crucial role in
the immune response. They
are primarily responsible
for producing and secreting
antibodies, also known as
immunoglobulins, which are
proteins that recognize and
bind to specific foreign
substances called antigens.
When the immune system
encounters an antigen, such
as a virus, bacteria, or
other foreign invader,
certain immune cells,
particularly B cells,
undergo a process called
activation. During
activation, B cells
differentiate into plasma
cells. This differentiation
process involves extensive
changes in gene expression
and cellular morphology,
transforming the B cell into
a highly specialized
antibody-producing factory.
Plasma cells are
characterized by their
abundant rough endoplasmic
reticulum (ER), which is
responsible for synthesizing
and processing large
quantities of antibodies.
These antibodies are
tailored to recognize and
bind to specific antigens
with high affinity. Once
produced, antibodies can
neutralize pathogens, mark
them for destruction by
other immune cells, or
facilitate their clearance
from the body through
various mechanisms.
Plasma cells are
particularly important
during the adaptive immune
response, which is
characterized by the
production of specific
antibodies tailored to the
invading pathogen. After the
initial encounter with an
antigen, plasma cells
continue to secrete
antibodies for a period of
time, providing the body
with a heightened defense
against reinfection. This
process is essential for
establishing immunity and
protecting against future
encounters with the same
pathogen.)");

const ScriptTemplate * ScriptTemplate::PlasmaCells() {return &PlasmaCellsScriptTemplate;}

constexpr ScriptTemplate PositiveGeneRegulationScriptTemplate("PositiveGeneRegulation.py", "\x01" R"(Some operons are also
subject to positive control
through a stimulatory
protein, such as catabolite
activator protein (CAP), an
activator of transcription.
When glucose (a preferred
food source of E. coli) is
scarce, CAP is activated by
binding with cyclic AMP
(cAMP). Activated CAP
attaches to the promoter of
the lac operon and increases
the affinity of RNA
polymerase, thus
accelerating transcription.
When glucose levels
increase, CAP detaches from
the lac operon, and
transcription returns to a
normal rate. CAP helps
regulate other operons that
encode enzymes used in
catabolic pathways)");

const ScriptTemplate * ScriptTemplate::PositiveGeneRegulation() {return &PositiveGeneRegulationScriptTemplate;}

constexpr ScriptTemplate ProliferationBTCellsScriptTemplate("ProliferationBTCells.py", "\x01" R"(aka spreading
In the body there are few
lymphocytes with antigen
receptors for any particular
epitope. In the lymph nodes,
an antigen is exposed to a
steady stream of lymphocytes
until a match is made. This
binding of a mature
lymphocyte to an antigen
initiates events that
activate the lymphocyte.
Once activated, a B or T
cell undergoes multiple cell
divisions. This
proliferation of lymphocytes
is called clonal selection.
Two types of clones are
produced: short lived
activated effector cells
that act immediately against
the antigen and long lived
memory cells that can give
rise to effector cells if
the same antigen is
encountered again.)");

const ScriptTemplate * ScriptTemplate::ProliferationBTCells() {return &ProliferationBTCellsScriptTemplate;}

constexpr ScriptTemplate ProtonMotiveForceScriptTemplate("ProtonMotiveForce.py", "\x01" R"(The proton motive force
(PMF) plays a crucial role
in the synthesis of ATP in
the mitochondria, a process
known as oxidative
phosphorylation due to the
electron transport chain. )");

const ScriptTemplate * ScriptTemplate::ProtonMotiveForce() {return &ProtonMotiveForceScriptTemplate;}

constexpr ScriptTemplate PublicHealthScriptTemplate("PublicHealth.py", "\x01" R"(Public health agencies form
a network, along with
private and voluntary
entities.
They aim to: prevent
epidemics and the spread of
disease, protect against
environmental hazards,
respond to disasters and
assists communities in
recovery, prevent injuries,
promote healthy behaviors,
assures the quality and
accessibility of health
services.
Definition by Sir Donald
Acheson: health of
communities, the science and
art of preventing disease,
prolonging life and
promoting, protecting and
improving health through the
organized efforts of
society.
Examples: vaccination,
fluoridation of drinking
water seatbelt policies,
non-smoking laws.
Investigated
problems/events:
environmental exposures
(lead and heavy metals, air
pollutants and other asthma
triggers), infectious
diseases (foodborne illness,
influenza and pneumonia),
injuries (increased
homicides, domestic
violence), non-infectious
diseases (increase in a
major birth defect), natural
disasters (hurricanes,
earthquakes), terrorism.
Work of epidemiologists: who
is sick? what are their
symptoms? when did they get
sick? where could they have
been exposed?
Mathematical models can
project how infectious
diseases progress to show
the likely outcome of an
epidemic(including in
plants) and help inform
public health and plant
health interventions.
Models use basic assumptions
or collected statistics
along with mathematics to
find parameters for various
infectious diseases and use
those parameters to
calculate the effects of
different interventions,
like mass vaccination
programs. The modelling can
help decide which
intervention(s) to avoid and
which to trial, or can
predict future growth
patterns, etc.)");

const ScriptTemplate * ScriptTemplate::PublicHealth() {return &PublicHealthScriptTemplate;}

constexpr ScriptTemplate RecognitionAndResponseScriptTemplate("RecognitionAndResponse.py", "\x01" R"(PATHOGENS - agents that
cause disease, infect a wide
range of animals, including
humans.
All animals have INNATE
IMMUNITY, a defense active
immediately upon infection.
Vertebrates also have
ADAPTIVE IMMUNITY.
Innate immunity (all
animals):
- Recognition of traits
shared by broad ranges of
pathogens, using a small set
of receptors.
- Rapid response
- Barrier defenses: skin,
mucous membranes, secretions

- Internal defenses:
phagocytic cells, natural
killer cells, antimicrobial
proteins, inflammatory
response.
Adaptive immunity
(vertebrates only):
- Recognition of traits
specific to particular
pathogens, using a vast
array of receptors.
- Slower response
- Humoral response:
antibodies defend against
infection in body fluids.
- Cell mediated response:
cytotoxic cells defend
against infection in body
cells.
Innate immunity is present
before any exposure to
pathogens and is effective
from the time of birth. it
involves nonspecific
responses to pathogens,
consists of external barrier
plus internal cellular and
chemical defenses.
Adaptive (acquired)
immunity, develops after
exposure to agent such as
microbes, toxins, or other
foreign substances. It
involves a very specific
response to pathogens.
)");

const ScriptTemplate * ScriptTemplate::RecognitionAndResponse() {return &RecognitionAndResponseScriptTemplate;}

constexpr ScriptTemplate RegulationOfTranscriptionInitiationScriptTemplate("RegulationOfTranscriptionInitiation.py", "\x01" R"(Chromatin-modifying enzymes
provide initial control of
gene expression by making a
region of DNA either more or
less able to bind the
transcription machinery.

Organizytion of a Typical
Eukaryotic Gene: Associated
with most eukaryotic genes
are multiple control
elements, segments of
noncoding DNA that serve as
binding sites for
transcription factors that
help regulate transcription.
Control elements and the
transcription factors they
bind are critical to the
precise regulation of gene
expression in different cell
types.

Roles of Transcription
Factors: To initiate
transcription, eukaryotic
RNA polymerase requires the
assistance of proteins
called transcription
factors. General
transcription factors are
essential for the
transcription of all
protein-coding genes. In
eukaryotes, high levels of
transcription of particular
genes depend on control
elements interacting with
specific transcription
factors.

Enhancers and Specific
Transcription Factors:
Proximal control elements
are located close to the
promoter. Distal control
elements, groupings of which
are called enhancers, may be
far away from a gene or even
located in an intron. An
activator is a protein that
binds to an enhancer and
stimulates transcription of
a gene. Activators have two
domains, one that binds DNA
and a second that activates
transcription. Bound
activators facilitate a
sequence of protein-protein
interactions that result in
transcription of a given
gene. Some transcription
factors function as
repressors, inhibiting
expression of a particular
gene by a variety of
methods. Some activators and
repressors act indirectly by
influencing chromatin
structure to promote or
silence transcription.)");

const ScriptTemplate * ScriptTemplate::RegulationOfTranscriptionInitiation() {return &RegulationOfTranscriptionInitiationScriptTemplate;}

constexpr ScriptTemplate SelfToleranceBTCellsScriptTemplate("SelfToleranceBTCells.py", "\x01" R"(Antigen receptors are
generated by random
rearrangement of DNA. As
lymphocytes mature in bone
marrow or the thymus, they
are tested for self
reactivity. Some B and T
cells with receptors
specific for the body’s own
molecules are destroyed by
apoptosis, or programmed
cell death. The remainder
are rendered nonfunctional.)");

const ScriptTemplate * ScriptTemplate::SelfToleranceBTCells() {return &SelfToleranceBTCellsScriptTemplate;}

constexpr ScriptTemplate TCellsAfterSelectionScriptTemplate("TCellsAfterSelection.py", "\x01" R"(Roles: T-cells are involved
in cell-mediated immunity,
where they directly interact
with infected or abnormal
cells. There are several
types of T-cells, including
helper T-cells, cytotoxic
T-cells, regulatory T-cells,
and memory T-cells. Helper
T-cells assist other immune
cells by secreting cytokines
and promoting their
activation. Cytotoxic
T-cells directly kill
infected or abnormal cells.
Regulatory T-cells suppress
immune responses to maintain
immune tolerance and prevent
autoimmune reactions. Memory
T-cells provide long-term
immunity upon re-exposure to
specific antigens.
Origin: T-cells also arise
from hematopoietic stem
cells in the bone marrow.
However, they undergo
further differentiation and
maturation in the thymus
gland. Immature T-cell
precursors migrate from the
bone marrow to the thymus,
where they undergo positive
and negative selection
processes that shape their
T-cell receptor (TCR)
repertoire and ensure
self-tolerance. After
selection, mature T-cells
migrate from the thymus to
secondary lymphoid organs,
where they encounter
antigens and become
activated.)");

const ScriptTemplate * ScriptTemplate::TCellsAfterSelection() {return &TCellsAfterSelectionScriptTemplate;}

constexpr ScriptTemplate Transcription1stScriptTemplate("Transcription1st.py", "\x01" R"(A GENE is a continuous
string of nucleotides
containing a region that
codes for an RNA molecule.
This region begins with a
PROMOTER and ends in
TERMINATOR.
For some genes, the encoded
RNA is used to to synthesize
a protein, in a process
called GENE EXPRESSION,
which could be divided in
two parts: TRANSCRIPTION and
TRANSLATION.
During transcription, the
DNA in the gene is used as a
template to make a messenger
RNA strand with the help of
the enzyme RNA POLYMERASE.
This process occurs in three
stages: INITIATION,
ELONGATION and TERMINATION.
Binding of RNA POLYMERASE
causes the DNA double helix
to unwind and open. During
ELONGATION RNA polymerase
slides along the template
DNA strand and links
nucleotides to the 5' end of
the growing RNA molecule.
When RNA polymerase reaches
TERMINATOR point DNA strand
and the messenger RNA
dissociate from each other.
Splicing:
The strand of messenger RNA
that is made during
transcription includes
regions called exons that
code for a protein, and
non-coding sections called
introns. Introns need to be
removed. This process is
called INTRON SPLICING and
is performed by a complex
made up of proteins and RNA
called a spliceosome. It
removes the intron segments
and joins the adjacent exons
to produce a mature
messenger RNA strand that
can leave the nucleus
through a nuclear pore and
enter the cytoplasm to begin
translation.
)");

const ScriptTemplate * ScriptTemplate::Transcription1st() {return &Transcription1stScriptTemplate;}

constexpr ScriptTemplate Translation2ndScriptTemplate("Translation2nd.py", "\x01" R"(The nitrogenous bases are
grouped into three letter
codes called CODONS.
Translation begins with the
messenger RNA strand binding
to the small ribosomal
subunit upstream of the
start codon. Each amino acid
is brought to the ribosome
by a specific transfer RNA
molecule (tRNA). The type of
amino acid is determined by
the anticodon sequence of
the transfer RNA.
After the tRNA binds to the
start codon, the large
ribosomal subunit binds to
form the translation complex
and INITIATION is complete. 
In the large ribosomal
subunit, there are three
distinct regions, called the
E, P, A sites. During
ELONGATION, individual amino
acids are brought to the
mRNA strand by a tRNA
through complementary base
pairing of the codons and
anticodons.
A charged tRNA molecule
binds to the A site and a
peptide bond forms.
The complex slides down
(towards 3') one codon to
the right where the now
uncharged tRNA exits from
the E site and the A site is
open to accept the next
transfer RNA. Elongation
will continue until stop
codon.
A release factor binds to
the A side at a stop codon
and the polypeptide is
released from the transfer
RNA in the P site. The
entire complex dissociates
and can reassemble.)");

const ScriptTemplate * ScriptTemplate::Translation2nd() {return &Translation2ndScriptTemplate;}

constexpr ScriptTemplate TryptophanOperonScriptTemplate("TryptophanOperon.py", "\x01" R"(The operon regulates and
expresses enzymes the
synthesise the amino acid
tryptophan in prokaryotes.
Repressor status when
synthesised by the regulator
gene: Tryptophan is not
usually available on the
extracellular region.
Prokaryotes need to
synthesise it. The repressor
is usually synthesised
INACTIVE
Corepressor/Inducer:
COREPRESSOR Trp = Tryptophan
Is the operon usually ON or
OFF?: ON
Status change What condition
triggers a status change?:
Trp is available
Corepressor/Inducer
role Where does the
corepressor/Inducer binds?
What is causes?: Trp
(corepressor) binds to the
allosteric site of the
repressor changing its
conformation into the ACTIVE
form.)");

const ScriptTemplate * ScriptTemplate::TryptophanOperon() {return &TryptophanOperonScriptTemplate;}

constexpr ScriptTemplate TypesOfChromosomeMutationsScriptTemplate("TypesOfChromosomeMutations.py", "\x01" R"(NEED TO WRITE: 
Point mutation located on __
triplet, __ base of the DNA
molecule, corresponding a
substitution/insertion/deletion.
The primary structure of the
protein is/is not affected
as the sequence of amino
acids is __ . Instead of __
there is __ . This causes a
samesense/missense/nonsense/frame
shift.
The final protein
conformation is
not/slightly/significantly
changed and this will/will
not affect the protein
function. (if the change
modifies the protein active
site)

DNA -> RNA -> Polypeptides

Deletion = A sequence of a
chromosome is missing.
Genetic material is lost.
ABC-DEF -> AC-DEF

Duplication = A sequence of
a chromosome is duplicated.
Genetic
material is gained
ABC-DEF -> ABBC-DEF

Inversion/reversal = A
sequence of a chromosome is
present in a wrong
sequence. There is no loss
or gain of genetic material.
ABC-DEF -> AED-CBF

Reciprocal/balanced
Translocation = 2
non-homologous chromosomes
exchange sequences. There is
no loss or gain of genetic
material.
ABC-DEF -> ABC-JKL
GH-IJKL -> GH-IDEF

Non-reciprocal translocation
or insertion - A chromosome
donates a sequence to
another non homologous
chromosome. There is no loss
or gain of genetic material.
ABC-DEF -> ABC-D
GH-IJKL -> GH-IEFJKL

Fusion - Two chromosomes
fuse together forming 1
chromosome. There is no loss
or gain of genetic material.
ABC-DEF -> ABC-DEFGH-IJKL
GH-IJKL -> 

Fission - A chromosomes
splits into 2 chromosomes.
There is no loss or
gain of genetic material.
ABC-DEF -> AB   C-DEF 

Substitution/insertion is
probably also possible

It can lead to loss or gain
of genetic material, or even
disruption of the gene
regulation.)");

const ScriptTemplate * ScriptTemplate::TypesOfChromosomeMutations() {return &TypesOfChromosomeMutationsScriptTemplate;}

constexpr ScriptTemplate TypesOfMutationsScriptTemplate("TypesOfMutations.py", "\x01" R"(Mutations produced by
germline cells, due to
meiotic nondisjunction.
Gene mutations (point
mutations)
Chromosome mutations (part
of a chromosome)
Genome mutations (whole
chromosome)
When are mutations produced?
Gene mutations and
chromosome mutations
->Happen during replication.
Genome mutations (whole
chromosome) -> Happen during
meiosis.
OR
Same sense: The mutated
triplet codes for the same
amino acid.
Mis-sense mutation: The
mutated triplet codes for a
different amino acid.
Non-sense mutation: The
mutated triplet codes for a
STOP codon.)");

const ScriptTemplate * ScriptTemplate::TypesOfMutations() {return &TypesOfMutationsScriptTemplate;}


// --------------------------------------------------------------------------------------------
constexpr ScriptTemplate emptyScriptTemplate(".py", "\x01" R"(from math import *
)");

const ScriptTemplate * ScriptTemplate::Empty() {
  return &emptyScriptTemplate;
}

}

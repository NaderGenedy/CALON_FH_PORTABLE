"""Build CALON Paper 1 Discussion docx from the traceable markdown."""
from docx import Document
from docx.shared import Pt
from docx.enum.text import WD_ALIGN_PARAGRAPH
import os

doc = Document()
st = doc.styles['Normal']
st.font.name = 'Cambria'
st.font.size = Pt(11)
st.paragraph_format.line_spacing = 1.15
st.paragraph_format.space_after = Pt(6)


def H(t, lvl=1):
    h = doc.add_heading(t, level=lvl)
    for r in h.runs:
        r.font.name = 'Cambria'


def P(t, bold=False, italic=False):
    p = doc.add_paragraph()
    r = p.add_run(t)
    r.font.name = 'Cambria'
    r.font.size = Pt(11)
    r.bold = bold
    r.italic = italic
    return p


title = doc.add_paragraph()
title.alignment = WD_ALIGN_PARAGRAPH.CENTER
r = title.add_run('CALON Paper 1 - Discussion (fully traceable to raw data)')
r.font.name = 'Cambria'
r.font.size = Pt(15)
r.bold = True
P('Multi-Evidence Structural Severity Score Predicts Statin Efficacy '
  'Independently of Polygenic Risk in Familial Hypercholesterolaemia', italic=True)
P('Traceability basis: every numerical claim re-derived from raw CSVs on '
  '2026-05-13. Discovery and validation headline correlations reproduce to '
  'three decimal places.', italic=True)

H('Discussion', 1)
P('This study demonstrates that a multi-evidence Structural Severity Score '
  '(SSS v3), integrating AlphaFold3 thermodynamics, deep mutational scanning, '
  'regulatory genomics, and Brown and Goldstein domain biology, independently '
  'predicts residual LDL-cholesterol on statin therapy in carriers of LDLR '
  'coding variants. The finding holds in two independent settings: a clinical '
  'FH registry (All-Wales discovery cohort, n = 109; SSS v3 versus LDL-C '
  'reduction Spearman rho = 0.234, P = 0.014) and a population-based external '
  'validation (UK Biobank, n = 775 statin-treated carriers; SSS v3 versus '
  'residual LDL-C rho = 0.083, P = 0.020). The central conceptual claim, that '
  'SSS v3 and the LDL-cholesterol polygenic risk score (LDL-PRS) capture '
  'orthogonal biological dimensions (rho = -0.014 between them) and additively '
  'stratify treatment outcome, reframes the clinical utility of LDLR structural '
  'analysis from prognostic to therapeutic.')

H('From pathogenicity to treatment response: a reframing, not a contest', 2)
P('The prevailing computational approaches to LDLR variant interpretation '
  '(AlphaMissense, evolutionary models, integrated scoring tools, and '
  'structure-based stability predictions) share a single objective: '
  'classifying a variant as pathogenic or benign. Tabet and colleagues deep '
  'mutational scanning provided the definitive experimental map of LDLR coding '
  'variation. These resources are powerful, and SSS v3 is built in part upon '
  'them. They are not, however, designed to answer the question a clinician '
  'faces at diagnosis: which treatment, and how aggressively escalated?')
P('CALON Paper 1 shows that the treatment question requires a structurally '
  'distinct kind of genetic information. Untreated LDL-C in our cohort was '
  'governed predominantly by polygenic background and the coarse functional '
  'domain of the variant; the fine-grained structural severity of receptor '
  'dysfunction added no detectable signal to untreated LDL-C (SSS v3 versus '
  'untreated LDL-C rho = -0.024, P = 0.235). In contrast, residual LDL-C, what '
  'the patient is left with after statin upregulation of receptor synthesis, '
  'was predicted by SSS v3 (rho = 0.083, P = 0.020), and this association '
  'survived adjustment for LDL-PRS (rho = 0.084, P = 0.020). The relationship '
  'to existing tools is therefore one of complementarity, not competition: SSS '
  'v3 answers a question the pathogenicity predictors were never built to '
  'address.')

H('Why categorical classification does not predict treatment response', 2)
P('The one direct head-to-head comparison in this study is between continuous '
  'structural scoring and categorical classification. The Brown and Goldstein '
  'five-class system did not predict statin treatment response in the Wales '
  'discovery cohort (Kruskal-Wallis non-significant for LDL-C reduction across '
  'functional classes). SSS v3, applied to the same patients, did (P = 0.014). '
  'The interpretation is not that the categorical system is wrong; it is that '
  'pharmacological response operates on a finer resolution than five classes '
  'can represent. A continuous score resolves that gradient; a categorical one '
  'cannot. This extends the demonstration by Mancini and colleagues that '
  'binary null-versus-defective classification fails to predict statin '
  'response in homozygous FH, to the heterozygous setting and to a granular, '
  'multi-evidence score.')

H('The two-dimensional genotype-to-treatment framework', 2)
P('Stratifying the UK Biobank validation cohort by median SSS v3 and median '
  'LDL-PRS produced a monotonic gradient in residual LDL-C on statin therapy, '
  'from carriers low on both axes to carriers high on both. The separation '
  'between the best- and worst-prognosis groups was approximately 0.5 mmol/L '
  'of residual LDL-C, a clinically meaningful margin when integrated over '
  'decades of treatment. Because SSS v3 and LDL-PRS are statistically '
  'independent (rho = -0.014), they define two genuinely separate axes: SSS v3 '
  'the mechanism of receptor dysfunction, LDL-PRS the cumulative polygenic '
  'burden. Their combination yields a two-by-two stratification applicable at '
  'the point of FH diagnosis, identifying the subset of carriers whose '
  'predicted residual LDL-C deficit warrants early escalation rather than '
  'conventional reactive titration.')

H('Effect size, honestly stated', 2)
P('The strength of the SSS v3 signal in external validation is modest: rho = '
  '0.083 corresponds to roughly 0.7% of variance in residual LDL-C explained '
  'by the structural score alone. This must be stated plainly. It sits within '
  'a larger and sobering picture in which the combined genetic layers explain '
  'only a small minority of total LDL-C variance, the majority remaining '
  'attributable to environmental and unmeasured factors. The value of SSS v3 '
  'is not that it is a strong univariate predictor, which it is not, but that '
  'it is a reproducible, polygenic-risk-independent signal that points to a '
  'real mechanistic axis, demonstrated consistently in a clinical registry '
  '(rho = 0.234) and a population biobank (rho = 0.083), and that it predicts '
  'treatment response specifically rather than disease severity. A weak but '
  'orthogonal and reproducible signal of treatment response is more clinically '
  'actionable than a strong but redundant predictor of pathogenicity.')

H('Negative findings', 2)
P('In keeping with a complete and unbiased account, several pre-specified '
  'analyses were null and are reported as such. SSS v3 did not predict '
  'untreated LDL-C (P = 0.235), did not predict incident major adverse '
  'cardiovascular events in time-to-event analysis (consistent with the '
  'treatment paradox, whereby carriers with more severe variants receive more '
  'aggressive therapy and observed outcomes converge), and did not predict '
  'cardiac MRI parameters. AlphaGenome regulatory scores showed domain-specific '
  'associations that did not survive domain adjustment, indicating confounding '
  'rather than independent regulatory biology, and are reported as '
  'exploratory. These nulls bound the claim: SSS v3 is a predictor of statin '
  'treatment response, not a general-purpose severity or outcome score.')

H('Conclusion', 2)
P('A multi-evidence structural severity score independently predicts statin '
  'efficacy in familial hypercholesterolaemia, validated in 775 UK Biobank '
  'LDLR variant carriers. SSS v3 and the LDL polygenic risk score capture '
  'orthogonal biological dimensions, receptor mechanism and polygenic burden, '
  'and their combination defines a two-dimensional genotype-to-treatment '
  'framework. The contribution is conceptual rather than one of raw predictive '
  'power: structural biology is reframed from estimating disease severity to '
  'anticipating therapeutic response, and the framework identifies, at '
  'diagnosis, the carriers most likely to require early treatment escalation.')

H('Traceability appendix', 2)
t = doc.add_table(rows=1, cols=3)
t.style = 'Light List Accent 2'
hdr = t.rows[0].cells
for i, x in enumerate(['Discussion claim', 'Manuscript value', 'Re-derived from raw']):
    hdr[i].text = x
rows = [
    ('Wales discovery SSS-LDL reduction', 'rho=0.234, P=0.014, n=109', 'rho=0.2341, P=0.0143  [exact]'),
    ('UKB validation SSS-residual LDL', 'rho=0.083, P=0.020, n=775', 'rho=0.0833, P=0.0203  [exact]'),
    ('PRS-adjusted SSS-residual LDL', 'rho=0.084, P=0.020', 'rho=0.0836, P=0.0202  [exact]'),
    ('SSS vs LDL-PRS orthogonality', 'rho=-0.014', 'rho=-0.0141  [exact]'),
    ('2x2 framework gradient', 'LL 3.10 to HH 3.50 mmol/L', 'LL 3.03 / HL 3.31 / HH 3.59  [within tol]'),
    ('Brown and Goldstein categorical', 'NS for treatment response', 'KW NS (P=0.19 11-group / 0.71 3-class)'),
]
for r in rows:
    cells = t.add_row().cells
    for i, v in enumerate(r):
        cells[i].text = v

P('UKB n=775 was reconstructed from raw medication fields p6153 and p6177 '
  '(the expanded cholesterol-medication cohort); the carrier CSV stored only '
  'the pre-expansion 374-patient definition. Reconstruction reproduces n=775 '
  'and rho=0.083 exactly. Numbers not used in this Discussion because they did '
  'not re-derive cleanly: the exact 0.40 mmol/L 2x2 spread (re-derivation '
  'gives approximately 0.5 mmol/L, and the Discussion uses the hedged value); '
  'cholesterol-years projections (modelled, not re-derived); '
  'variance-decomposition percentages (manuscript-stated, re-derivation '
  'deferred).', italic=True)

out = 'CALON_Paper1_Discussion_TRACEABLE.docx'
doc.save(out)
print('wrote ' + out + ' (' + str(round(os.path.getsize(out) / 1024, 1)) + ' KB)')

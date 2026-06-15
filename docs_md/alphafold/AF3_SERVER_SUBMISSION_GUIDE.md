# AlphaFold3 Server Submission Guide

**Date:** 2026-03-17
**Project:** CALON-Structure Phase 2a
**Server:** https://alphafoldserver.com
**Jobs required:** 5 (out of your 30/day limit)

---

## Job 1: LDLR Wild-Type (Monomer)

- **Job name:** `CALON_LDLR_wildtype`
- **Entity type:** Protein
- **UniProt:** P01130
- **Copies:** 1

**Sequence (860 residues):**
```
MGPWGWKLRWTVALLLAAAGTAVGDRCERNEFQCQDGKCISYKWVCDGSAECQDGSDESQETCLSVTCKSGDFSCGGRVNRCIPQFWRCDGQVDCDNGSDEQGCPPKTCSQDEFRCHDGKCISRQFVCDSDRDCLDGSDEASCPVLTCGPASFQCNSSTCIPQLWACDNDPDCEDGSDEWPQRCRGLYVFQGDSSPCSAFEFHCLSGECIHSSWRCDGGPDCKDKSDEENCAVATCRPDEFQCSDGNCIHGSRQCDREYDCKDMSDEVGCVNVTLCEGPNKFKCHSGECITLDKVCNMARDCRDWSDEPIKECGTNECLDNNGGCSHVCNDLKIGYECLCPDGFQLVAQRRCEDIDECQDPDTCSQLCVNLEGGYKCQCEEGFQLDPHTKACKAVGSIAYLFFTNRHEVRKMTLDRSEYTSLIPNLRNVVALDTEVASNRIYWSDLSQRMICSTQLDRAHGVSSYDTVISRDIQAPDGLAVDWIHSNIYWTDSVLGTVSVADTKGVKRKTLFRENGSKPRAIVVDPVHGFMYWTDWGTPAKIKKGGLNGVDIYSLVTENIQWPNGITLDLLSGRLYWVDSKLHSISSIDVNGGNRKTILEDEKRLAHPFSLAVFEDKVFWTDIINEAIFSANRLTGSDVNLLAENLLSPEDMVLFHNLTQPRGVNWCERTTLSNGGCQYLCLPAPQINPHSPKFTCACPDGMLLARDMRSCLTEAEAAVATQETSTVRLKVSSTAVRTQHTTTRPVPDTSRLPGATPGLTTVEIVTMSHQALGDVAGRGNEKKPSSVRALSIVLPIVLLVFLCLGVFLLWKNWRLKNINSINFDNPVYQKTTEDEVHICHNQDGYSYPSRQMVSLEDDVA
```

**Purpose:** Primary structure for ALL LDLR variant modelling via FoldX. This is the most important job.

---

## Job 2: PCSK9 Wild-Type (Monomer)

- **Job name:** `CALON_PCSK9_wildtype`
- **Entity type:** Protein
- **UniProt:** Q8NBP7
- **Copies:** 1

**Sequence (692 residues):**
```
MGTVSSRRSWWPLPLLLLLLLLLGPAGARAQEDEDGDYEELVLALRSEEDGLAEAPEHGTTATFHRCAKDPWRLPGTYVVVLKEETHLSQSERTARRLQAQAARRGYLTKILHVFHGLLPGFLVKMSGDLLELALKLPHVDYIEEDSSVFAQSIPWNLERITPPRYRADEYQPPDGGSLVEVYLLDTSIQSDHREIEGRVMVTDFENVPEEDGTRFHRQASKCDSHGTHLAGVVSGRDAGVAKGASMRSLRVLNCQGKGTVSGTLIGLEFIRKSQLVQPVGPLVVLLPLAGGYSRVLNAACQRLARAGVVLVTAAGNFRDDACLYSPASAPEVITVGATNAQDQPVTLGTLGTNFGRCVDLFAPGEDIIGASSDCSTCFVSQSGTSQAAAHVAGIAAMMLSAEPELTLAELRQRLIHFSAKDVINEAWFPEDQRVLTPNLVAALPPSTHGAGWQLFCRTVWSAHSGPTRMATAVARCAPDEELLSCSSFSRSGKRRGERMEAQGGKLVCRAHNAFGGEGVYAIARCCLLPQANCSVHTAPPAEASMGTRVHCHQQGHVLTGCSSHWEVEDLGTHKPPVLRPRGQPNQCVGHREASIHASCHCHAPGLECKVKEHGIPAPQEQVTVACEEGWTLTGCSALPGTSHVLGAYAVDNTCVVRSRDVSTTGSTSEGAVTAVAICCRSRHLAQASQELQ
```

**Purpose:** Structure for PCSK9 gain-of-function variant modelling.

---

## Job 3: ApoB-100 Receptor-Binding Domain (Monomer)

- **Job name:** `CALON_ApoB_RBD`
- **Entity type:** Protein
- **UniProt:** P04114 (residues 3359-3600 only)
- **Copies:** 1

**Sequence (242 residues):**
```
QSDIVAHLLSSSSSVIDALQYKLEGTTRLTRKRGLKLATALSLSNKFVEGSHNSTVSLTTKNMEVSVATTTKAQIPILRMNFKQELNGNTKSKPTVSSSMEFKYDFNSSMLYSTAKGAVDHKLSLESLTSYFSIESSTKGDVKGSVLSREYSGTIASEANTYLNSKSTRSSVKLQGTSKIDDIWNLEVKENFAGEATLQRIYSLWEHSTKNHLQLEGLFFTNGEHTSKATLELSPWQMSALV
```

**Purpose:** Structure for APOB:c.10580G>A (p.R3527Q) and other APOB RBD variants. R3527 falls at position 169 within this fragment (3527 - 3359 + 1 = 169).

**Note:** Full ApoB-100 is 4,563 residues — too large for AF3. The RBD is sufficient because all known FH-causing APOB variants cluster in this region.

---

## Job 4: LDLR + PCSK9 Complex (Multimer)

- **Job name:** `CALON_LDLR_PCSK9_complex`
- **Entity type:** Protein (2 chains)

**Chain A — LDLR (use same sequence as Job 1):**
- Copies: 1

**Chain B — PCSK9 (use same sequence as Job 2):**
- Copies: 1

**Purpose:** Model the LDLR-PCSK9 degradation pathway interface. Variants near this interface may enhance PCSK9-mediated LDLR degradation (gain-of-function mechanism). Interface residues will be extracted for the "interface proximity" SSS metric.

---

## Job 5: LDLR + ApoB RBD Complex (Multimer)

- **Job name:** `CALON_LDLR_ApoB_complex`
- **Entity type:** Protein (2 chains)

**Chain A — LDLR (use same sequence as Job 1):**
- Copies: 1

**Chain B — ApoB RBD (use same sequence as Job 3):**
- Copies: 1

**Purpose:** Model the LDL particle binding interface. Variants disrupting this interface directly impair cholesterol clearance. Critical for interpreting both LDLR ligand-binding domain variants AND APOB:p.R3527Q.

---

## How to Submit on alphafoldserver.com

### For monomer jobs (Jobs 1-3):
1. Go to https://alphafoldserver.com
2. Click "New Job"
3. Set job name (e.g., `CALON_LDLR_wildtype`)
4. Add entity → Protein → Paste sequence
5. Set copies = 1
6. Click "Submit"

### For multimer jobs (Jobs 4-5):
1. Go to https://alphafoldserver.com
2. Click "New Job"
3. Set job name (e.g., `CALON_LDLR_PCSK9_complex`)
4. Add entity → Protein → Paste Chain A sequence → Copies = 1
5. Add another entity → Protein → Paste Chain B sequence → Copies = 1
6. Click "Submit"

### Important notes:
- AF3 Server outputs mmCIF format — you will need to convert to PDB for FoldX
- Download ALL output files (ranked models, confidence scores, PAE plots)
- The pLDDT scores in the output are needed for the SSS pLDDT component
- Multimer jobs may take longer than monomer jobs
- Save the JSON input files from the server for reproducibility

---

## After AF3 Jobs Complete

1. Download all 5 result sets (mmCIF + confidence files)
2. Save to `alphafold/structures/` directory
3. Convert mmCIF → PDB using `gemmi convert` or BioPython
4. Run FoldX RepairPDB on wild-type structures
5. Run FoldX BuildModel for each of the top 50 missense variants
6. Extract interface residues from complex structures (Jobs 4-5)

---

## Expected Timeline

- **Submission:** ~15 minutes (all 5 jobs)
- **Monomer results:** 1-4 hours typically
- **Multimer results:** 2-8 hours typically
- **All results available:** Same day (within your 30-job daily limit)

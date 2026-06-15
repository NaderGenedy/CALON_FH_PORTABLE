"""Extract numerical claims from a manuscript .docx using regex."""
import argparse, json, re, sys
from pathlib import Path

PATTERNS = {
    'AUC_with_CI':    r'AUC\s*[=:]?\s*([01]\.\d{2,4})\s*(?:\([^)]*?(\d\.\d{2,4})\s*[-\u2013]\s*(\d\.\d{2,4})[^)]*\))?',
    'NRI':            r'NRI\s*(?:of\s+\w+\s+over\s+\w+\s+was)?\s*[=:]?\s*([+-]?\d\.\d{2,4})',
    'IDI':            r'IDI\s*[=:]?\s*([+-]?\d\.\d{2,4})',
    'HR_with_CI':     r'HR\s*(?:[=:]\s*)?(\d+\.\d{2,3})\s*\(([\d.]+)\s*[-\u2013]\s*([\d.]+)\)',
    'OR_with_CI':     r'OR\s*(?:[=:]\s*)?(\d+\.\d{2,3})\s*\(([\d.]+)\s*[-\u2013]\s*([\d.]+)\)',
    'p_value':        r'p\s*[=<>]\s*(\d\.\d+(?:\s*[x\u00d7]\s*10\^?[-\u2212]?\d+)?)',
    'sample_size':    r'n\s*=\s*([\d,]+)',
    'percentage':     r'(\d+\.\d+)%',
    'brier':          r'Brier\s*(?:score)?\s*[=:]?\s*(\d\.\d{3,4})',
    'calib_slope':    r'[Cc]alibration\s+slope\s*(?:of)?\s*[=:]?\s*(\d+\.\d{2})',
    'sensitivity':    r'[Ss]ensitivity\s*[=:]?\s*(\d+\.\d+)%?',
    'specificity':    r'[Ss]pecificity\s*[=:]?\s*(\d+\.\d+)%?',
}

def load_text(path):
    p = Path(path)
    if p.suffix.lower() == '.docx':
        from docx import Document
        doc = Document(p)
        return '\n\n'.join(para.text for para in doc.paragraphs if para.text.strip())
    return p.read_text(encoding='utf-8', errors='replace')

def extract_claims(text):
    claims = []; seen = set()
    for kind, pat in PATTERNS.items():
        for m in re.finditer(pat, text, re.IGNORECASE):
            start = max(0, m.start() - 60); end = min(len(text), m.end() + 60)
            ctx = text[start:end].replace('\n', ' ').strip()
            key = (kind, m.group(0), ctx[:50])
            if key in seen: continue
            seen.add(key)
            try: val = float(m.group(1).replace(',', ''))
            except (ValueError, AttributeError): continue
            claims.append({'kind': kind, 'value': val, 'raw_match': m.group(0), 'context': ctx, 'position': m.start()})
    return claims

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('manuscript'); ap.add_argument('--out', default=None)
    args = ap.parse_args()
    text = load_text(args.manuscript); claims = extract_claims(text)
    out = {'source': str(args.manuscript), 'n_claims': len(claims), 'claims': claims}
    payload = json.dumps(out, indent=2, ensure_ascii=False)
    if args.out: Path(args.out).write_text(payload, encoding='utf-8'); print(f'Wrote {len(claims)} claims to {args.out}')
    else: print(payload)

if __name__ == '__main__': main()

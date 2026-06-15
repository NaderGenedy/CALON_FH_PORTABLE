@echo off
set PYTHONPATH=
set PYTHONHOME=
"C:\Program Files\ChimeraX 1.11.1\bin\ChimeraX-console.exe" --nogui --offscreen --exit --cmd "open 'C:\Users\nader\Downloads\calon_ukb_pipeline\alphafold\foldx\LDLR_wt_Repair.pdb' ; set bgColor white ; save 'C:\Users\nader\Downloads\calon_ukb_pipeline\alphafold\analysis\figures\chimerax\test_render.png' width 1000 height 800 supersample 2"

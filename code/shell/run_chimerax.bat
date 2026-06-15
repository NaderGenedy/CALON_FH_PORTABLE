@echo off
set PYTHONPATH=
set PYTHONHOME=
set CHIMERAX="C:\Program Files\ChimeraX 1.11.1\bin\ChimeraX-console.exe"
set FIGDIR=C:\Users\nader\Downloads\calon_ukb_pipeline\alphafold\analysis\figures\chimerax

echo === Running Figure CX01: Risk Atlas ===
%CHIMERAX% --nogui --offscreen --exit --cmd "open %FIGDIR%\figure_CX01.cxc"

echo === Running Figure CX02: pLDDT Map ===
%CHIMERAX% --nogui --offscreen --exit --cmd "open %FIGDIR%\figure_CX02.cxc"

echo === Running Figure CX03: PCSK9 Interface ===
%CHIMERAX% --nogui --offscreen --exit --cmd "open %FIGDIR%\figure_CX03.cxc"

echo === Running Figure CX04: ApoB Interface ===
%CHIMERAX% --nogui --offscreen --exit --cmd "open %FIGDIR%\figure_CX04.cxc"

echo === Running Figure CX05: Zero-Overlap ===
%CHIMERAX% --nogui --offscreen --exit --cmd "open %FIGDIR%\figure_CX05.cxc"

echo === Running Figure CX06: Domain Architecture ===
%CHIMERAX% --nogui --offscreen --exit --cmd "open %FIGDIR%\figure_CX06.cxc"

echo === Running Figure CX07: Critical Residues ===
%CHIMERAX% --nogui --offscreen --exit --cmd "open %FIGDIR%\figure_CX07.cxc"

echo === Running Figure CX08: Hotspot Heatmap ===
%CHIMERAX% --nogui --offscreen --exit --cmd "open %FIGDIR%\figure_CX08.cxc"

echo === Running Figure CX09: ddG Map ===
%CHIMERAX% --nogui --offscreen --exit --cmd "open %FIGDIR%\figure_CX09.cxc"

echo === Running Figure CX10: Drug Targets ===
%CHIMERAX% --nogui --offscreen --exit --cmd "open %FIGDIR%\figure_CX10.cxc"

echo === ALL DONE ===
dir /b "%FIGDIR%\Figure_CX*.png"

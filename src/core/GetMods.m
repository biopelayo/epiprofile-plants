function Mods = GetMods()
%% GetMods — Histone/PTM modification catalog (monoisotopic & average deltas)
%
% PLAIN-ENGLISH PURPOSE (non-experts)
%   Returns a structure 'Mods' that acts as a *lookup table* for post-translational
%   modifications (PTMs) and labeling reagents. Each entry defines:
%     • Mods.set{n}  : which residue letters are allowed for that modification
%                      (e.g., 'K', 'R', 'STY', etc.; 1-letter, uppercase).
%     • Mods.name{n} : short code for the modification (e.g., 'ac', 'me1', 'ph').
%     • Mods.mass(:,n): [mono; avge] mass *deltas* (Da) to be applied when the mod
%                       is present (monoisotopic in row 1, average in row 2).
%
% IMPORTANT
%   • Masses are *deltas relative to the unmodified residue in this workflow*.
%   • Some deltas are **derivatization-aware** (common in histone proteomics):
%     e.g., in propionylation workflows, K-me1 can be encoded as (propionyl + methyl).
%   • Isotopic labels (SILAC Lys8/Arg10, CD3/13CD3 methyls) are provided as separate
%     entries or combinatorial variants.
%
% FIELDS
%   Mods.mass : [2 × N] double  — row 1 = monoisotopic, row 2 = average
%   Mods.set  : 1 × N cellstr   — residue letters (where the mod is allowed)
%   Mods.name : 1 × N cellstr   — short name / code for the modification
%
% REFERENCES
%   Names inspired by common usage; masses compatible with Unimod-style deltas.
%   See: http://www.unimod.org/  (Cross-check when extending this table.)
%
% CAVEATS
%   • [Inference] Some heavy-label variants combine derivatization (+pr/+ac) and
%     methyl isotopologues (CD3 or 13CD3). Shortcodes like 'hme1', 'me11' reflect
%     different isotope patterns; see the README for decoding.
%   • This table is project-specific; adapt or extend with care and document changes.

    bmono = 1;              % row index for monoisotopic mass
    bavge = 2;              % row index for average mass
    nMod  = 500;            % rough preallocation capacity (number of entries)
    n     = 0;              % actual number of filled entries

    % Pre-allocate mass matrix (two rows: mono/avge; columns filled progressively)
    Mods.mass = zeros([2,nMod]);

    %% --- Core acylations / methylations (histone-centric) ---
    % pr — Propionylation (often derivatization of K/S/T/Y in histone workflows)
    n = n + 1;
    Mods.set{n}      = 'KSTY';
    Mods.name{n}     = 'pr';
    Mods.mass(bmono,n) = 56.026215;   % +C3H4O (propionyl) delta
    Mods.mass(bavge,n) = 56.0633;

    % ac — Acetylation
    n = n + 1;
    Mods.set{n}      = 'KSTY';
    Mods.name{n}     = 'ac';
    Mods.mass(bmono,n) = 42.010565;
    Mods.mass(bavge,n) = 42.0367;

    % me1 (K) — [Inference] K-monomethyl delta *in a propionyl workflow*:
    %           equals propionyl (+56.026215) + methyl (+14.01565) = 70.041865
    n = n + 1;
    Mods.set{n}      = 'K';
    Mods.name{n}     = 'me1';
    Mods.mass(bmono,n) = 70.041865;
    Mods.mass(bavge,n) = 70.0898;

    % Me1 (R) — Arginine monomethylation (note uppercase 'M' to distinguish legacy code)
    n = n + 1;
    Mods.set{n}      = 'R';
    Mods.name{n}     = 'Me1';
    Mods.mass(bmono,n) = 14.015650;   % classic +CH3 on arginine
    Mods.mass(bavge,n) = 14.0266;

    % me2 — Di-methylation (applicable to K or R as per set)
    n = n + 1;
    Mods.set{n}      = 'KR';
    Mods.name{n}     = 'me2';
    Mods.mass(bmono,n) = 28.031300;
    Mods.mass(bavge,n) = 28.0532;

    % me3 (K) — Tri-methylation (lysine)
    n = n + 1;
    Mods.set{n}      = 'K';
    Mods.name{n}     = 'me3';
    Mods.mass(bmono,n) = 42.046950;
    Mods.mass(bavge,n) = 42.0797;

    % ph — Phosphorylation (S/T/Y)
    n = n + 1;
    Mods.set{n}      = 'STY';
    Mods.name{n}     = 'ph';
    Mods.mass(bmono,n) = 79.966331;
    Mods.mass(bavge,n) = 79.9799;

    % cr — Crotonylation (K)
    n = n + 1;
    Mods.set{n}      = 'K';
    Mods.name{n}     = 'cr';
    Mods.mass(bmono,n) = 68.026215;
    Mods.mass(bavge,n) = 68.0740;

    % ub — Diglycine remnant on K (ubiquitylation signature after trypsin) +114.042927
    n = n + 1;
    Mods.set{n}      = 'K';
    Mods.name{n}     = 'ub';
    Mods.mass(bmono,n) = 114.042927;
    Mods.mass(bavge,n) = 114.1026;

    % ox — Oxidation (M → M[+16], e.g., Met sulfoxide)
    n = n + 1;
    Mods.set{n}      = 'M';
    Mods.name{n}     = 'ox';
    Mods.mass(bmono,n) = 15.994915;
    Mods.mass(bavge,n) = 15.9994;

    %% --- SILAC / heavy labels (Lys8, Arg10) ---
    % lak — Lys8 (SILAC) on K: +8.014199 (13C6 15N2-Lys vs light)
    n = n + 1;
    Mods.set{n}      = 'K';
    Mods.name{n}     = 'lak';
    Mods.mass(bmono,n) = 8.014199;
    Mods.mass(bavge,n) = 7.9427;

    % lar — Arg10 (SILAC) on R: +10.008269 (13C6 15N4-Arg vs light)
    n = n + 1;
    Mods.set{n}      = 'R';
    Mods.name{n}     = 'lar';
    Mods.mass(bmono,n) = 10.008269;
    Mods.mass(bavge,n) = 9.9296;

    %% --- Heavy acetyl / heavy methyl families (derivatization & isotopologues) ---
    % hac — [Inference] heavy acetyl on K (e.g., 13C2-acetyl) ~ +44.0173
    n = n + 1;
    Mods.set{n}      = 'K';
    Mods.name{n}     = 'hac';
    Mods.mass(bmono,n) = 44.017274;
    Mods.mass(bavge,n) = 44.0220;

    % hme1 — [Inference] K-me1 with CD3-methyl instead of CH3, in propionyl workflow
    %         me1(+70.041865) + (CD3−CH3 ≈ +3.01883) = +73.060695
    n = n + 1;
    Mods.set{n}      = 'K';
    Mods.name{n}     = 'hme1';
    Mods.mass(bmono,n) = 73.060695;
    Mods.mass(bavge,n) = 73.1083;

    % hme2 — two CD3 methyls vs CH3: me2(+28.0313) + 2×3.01883 = +34.06896
    n = n + 1;
    Mods.set{n}      = 'K';
    Mods.name{n}     = 'hme2';
    Mods.mass(bmono,n) = 34.06896;
    Mods.mass(bavge,n) = 34.0902;

    % hme3 — three CD3 methyls vs CH3: me3(+42.04695) + 3×3.01883 = +51.10344
    n = n + 1;
    Mods.set{n}      = 'K';
    Mods.name{n}     = 'hme3';
    Mods.mass(bmono,n) = 51.10344;
    Mods.mass(bavge,n) = 51.1352;

    %% --- Combos with SILAC Lys8 (lak) on K ---
    % prlak — propionyl + Lys8
    n = n + 1;
    Mods.set{n}      = 'K';
    Mods.name{n}     = 'prlak';
    Mods.mass(bmono,n) = 64.040414;   % 56.026215 + 8.014199
    Mods.mass(bavge,n) = 64.006;

    % aclak — acetyl + Lys8
    n = n + 1;
    Mods.set{n}      = 'K';
    Mods.name{n}     = 'aclak';
    Mods.mass(bmono,n) = 50.024764;   % 42.010565 + 8.014199
    Mods.mass(bavge,n) = 49.9794;

    % me1lak — [Inference] (propionyl + me1) + Lys8
    n = n + 1;
    Mods.set{n}      = 'K';
    Mods.name{n}     = 'me1lak';
    Mods.mass(bmono,n) = 78.056064;   % 70.041865 + 8.014199
    Mods.mass(bavge,n) = 78.0325;

    % me2lak — me2 + Lys8
    n = n + 1;
    Mods.set{n}      = 'K';
    Mods.name{n}     = 'me2lak';
    Mods.mass(bmono,n) = 36.045499;   % 28.031300 + 8.014199
    Mods.mass(bavge,n) = 35.9959;

    % me3lak — me3 + Lys8
    n = n + 1;
    Mods.set{n}      = 'K';
    Mods.name{n}     = 'me3lak';
    Mods.mass(bmono,n) = 50.061149;   % 42.046950 + 8.014199
    Mods.mass(bavge,n) = 50.0224;

    %% --- Isotopologue ladder for methyl groups using 13CD3 (vs CH3) on K ---
    % me11 — me1 with ONE 13CD3 methyl instead of CH3: +4.022185 over me1
    n = n + 1;
    Mods.set{n}      = 'K';
    Mods.name{n}     = 'me11';
    Mods.mass(bmono,n) = 74.06405;    % 70.041865 + 4.022185
    Mods.mass(bavge,n) = 74.1009;

    % me21 — me2 with ONE 13CD3 (and one CH3): +4.022185 over me2
    n = n + 1;
    Mods.set{n}      = 'K';
    Mods.name{n}     = 'me21';
    Mods.mass(bmono,n) = 32.053485;   % 28.031300 + 4.022185
    Mods.mass(bavge,n) = 32.0643;

    % me22 — me2 with TWO 13CD3: +2 × 4.022185 over me2
    n = n + 1;
    Mods.set{n}      = 'K';
    Mods.name{n}     = 'me22';
    Mods.mass(bmono,n) = 36.07567;    % 28.031300 + 8.04437
    Mods.mass(bavge,n) = 36.0754;

    % me31 — me3 with ONE 13CD3 (2×CH3 + 1×13CD3)
    n = n + 1;
    Mods.set{n}      = 'K';
    Mods.name{n}     = 'me31';
    Mods.mass(bmono,n) = 46.069135;   % 42.046950 + 4.022185
    Mods.mass(bavge,n) = 46.0908;

    % me32 — me3 with TWO 13CD3 (1×CH3 + 2×13CD3)
    n = n + 1;
    Mods.set{n}      = 'K';
    Mods.name{n}     = 'me32';
    Mods.mass(bmono,n) = 50.09132;    % 42.046950 + 8.04437
    Mods.mass(bavge,n) = 50.1019;

    % me33 — me3 with THREE 13CD3 (fully heavy methyls)
    n = n + 1;
    Mods.set{n}      = 'K';
    Mods.name{n}     = 'me33';
    Mods.mass(bmono,n) = 54.113505;   % 42.046950 + 12.066555
    Mods.mass(bavge,n) = 54.113;

    %% --- Building block for 13CD3 (legacy helper) ---
    % cd3 — (used as a primitive delta in some generators; here associated to 'M' for legacy)
    n = n + 1;
    Mods.set{n}      = 'M';
    Mods.name{n}     = 'cd3';
    Mods.mass(bmono,n) = 4.022185;    % delta of 13CD3 vs CH3 (legacy convenience)
    Mods.mass(bavge,n) = 4.0111;

    %% Finalize: shrink mass matrix to actual number of entries
    if n < nMod
        Mods.mass = Mods.mass(1:2,1:n);
    end
end

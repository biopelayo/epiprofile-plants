function His = init_histone0(special)
%% init_histone0
% PURPOSE (English, non-proteomics phrasing):
%   Build a curated catalog (struct "His") of histone peptides and their
%   fixed derivatization/PTM settings for targeted extraction/quantification.
%   Each peptide row includes: a human-readable name, the amino-acid sequence,
%   fixed modifications by position, a target charge state, and a computed
%   theoretical m/z.
%
% WHY THIS LAYOUT:
%   - Focus on H3/H4 anchors and biologically relevant windows (e.g., H3 27–40),
%     including Arabidopsis-specific variants (e.g., H3.3-like "TT" motif).
%   - Include optional PTM probes (e.g., S10ph, K9ac) to confirm known marks.
%   - Keep a minimal H2A block mostly as QC anchors (can be dropped if unused).
%
% FIELD SEMANTICS:
%   His.out_filename{no,1} : human label / output prefix (may be used for folders).
%   His.pep_seq{no,1}      : peptide sequence (string).
%   His.mod_type{no,1}     : fixed mods "pos,type;" grammar (0 = N-term; 1..N = residues).
%                            types: pr=propionylation, ph=phosphorylation, ac=acetylation.
%   His.pep_ch(no,1)       : target charge for m/z calculation (1/2/3).
%   His.pep_mz(no,1)       : monoisotopic m/z computed from sequence + mods + charge.
%   His.seq_godel(no,1)    : lightweight numeric fingerprint of "seq+mods" (optional).
%
% CAUTION:
%   Some "out_filename" labels are reused across variants (same prefix). If those labels
%   become folder names downstream, you may overwrite files. Consider suffixes if needed.

no = 0;

% checklist
no = no + 1;
His.out_filename{no,1} = 'H3_01_3_8';          % H3 short N-term anchor (residues 3–8)
His.pep_seq{no,1} = 'TKQTAR';                  % T K Q T A R
His.mod_type{no,1} = '0,pr;2,pr;';             % fixed propionylation at N-term (0) and K at pos 2
His.pep_ch(no,1) = 1;                           % short peptide: +1 charge is reasonable
His.pep_mz(no,1) = calculate_pepmz0(His,no,special); % compute theoretical m/z with current AA masses and mods
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];    % concat to form a unique signature string
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq))); % small numeric fingerprint (optional)

no = no + 1;
His.out_filename{no,1} = 'H3_02_9_17';         % H3 window 9–17: classic N-term tryptic segment
His.pep_seq{no,1} = 'KSTGGKAPR';               % includes K sites relevant around K9/K14 region
His.mod_type{no,1} = '0,pr;1,pr;6,pr;';        % propionylation at N-term and Lys positions
His.pep_ch(no,1) = 2;                           % typical +2 for this length
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_03_18_26';        % H3 window 18–26
His.pep_seq{no,1} = 'KQLATKAAR';
His.mod_type{no,1} = '0,pr;1,pr;6,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_04_27_40';        % H3 K27/K36 region (key biological hotspot)
His.pep_seq{no,1} = 'KSAPATGGVKKPHR';         % "AT" motif version (H3.1-like)
His.mod_type{no,1} = '0,pr;1,pr;10,pr;11,pr;'; % derivatized N-term + Lys positions
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H4_01_4_17';         % H4 tail 4–17: Lys-rich anchor (propionylation-friendly)
His.pep_seq{no,1} = 'GKGGKGLGKGGAKR';
His.mod_type{no,1} = '0,pr;2,pr;5,pr;9,pr;13,pr;'; % derivatization at N-term and multiple Lys
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H4_02_20_23';        % H4 short anchor 20–23
His.pep_seq{no,1} = 'KVLR';
His.mod_type{no,1} = '0,pr;1,pr;';             % N-term + Lys derivatization
His.pep_ch(no,1) = 1;                           % very short: +1
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

%------------------H3------------------
no = no + 1;
His.out_filename{no,1} = 'H3_02a_9_17';        % H3 9–17 with explicit phosphorylation probe
His.pep_seq{no,1} = 'KSTGGKAPR';
His.mod_type{no,1} = '0,pr;1,pr;2,ph;6,pr;';   % add S phosphorylation (pos 2) as a targeted check
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_02b_9_17';        % H3 9–17 with explicit acetylation probe
His.pep_seq{no,1} = 'KSTGGKAPR';
His.mod_type{no,1} = '0,pr;1,pr;2,ac;6,pr;';   % add K/S acetylation probe (as defined by get_mod_mass grammar)
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_05_41_49';        % H3 internal stable anchor (away from dense PTM regions)
His.pep_seq{no,1} = 'YRPGTVALR';
His.mod_type{no,1} = '0,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_06_53_63';        % H3 window 53–63 (plant-tuned start)
His.pep_seq{no,1} = 'KYQKSTELLIR';
His.mod_type{no,1} = '0,pr;1,pr;4,pr;';        % propionylation at N-term and selected Lys
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_06a_53_63';       % same window with leading R variant
His.pep_seq{no,1} = 'RYQKSTELLIR';
His.mod_type{no,1} = '0,pr;4,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_07_73_83';        % H3 window 73–83 (internal anchor)
His.pep_seq{no,1} = 'EIAQDFKTDLR';
His.mod_type{no,1} = '0,pr;7,pr;';             % N-term + Lys
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_08_117_128';      % H3 C-terminal anchor 117–128
His.pep_seq{no,1} = 'VTIMPKDIQLAR';
His.mod_type{no,1} = '0,pr;6,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_09u_64_135';      % group of micro-anchors within 64–135 (same label; beware of file collisions)
His.pep_seq{no,1} = 'KLPFQR';
His.mod_type{no,1} = '0,pr;1,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_09u_64_135';      % same label (second micro-anchor)
His.pep_seq{no,1} = 'RIRGER';
His.mod_type{no,1} = '0,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_09u_64_135';      % same label (third micro-anchor)
His.pep_seq{no,1} = 'IRGERA';
His.mod_type{no,1} = '0,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_11_3_8';          % N-term polymorphism set (same label reused)
His.pep_seq{no,1} = 'TKQTAR';
His.mod_type{no,1} = '0,pr;2,pr;';
His.pep_ch(no,1) = 1;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_11_3_8';          % same label (variant)
His.pep_seq{no,1} = 'TKQSAR';
His.mod_type{no,1} = '0,pr;2,pr;';
His.pep_ch(no,1) = 1;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_11_3_8';          % same label (variant)
His.pep_seq{no,1} = 'SNQTAR';
His.mod_type{no,1} = '0,pr;';
His.pep_ch(no,1) = 1;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_12_9_17';         % 9–17 variant set (label reused across variants)
His.pep_seq{no,1} = 'KSTGGKAPR';
His.mod_type{no,1} = '0,pr;1,pr;6,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_12_9_17';         % same label (variant)
His.pep_seq{no,1} = 'KSTGGKGPR';
His.mod_type{no,1} = '0,pr;1,pr;6,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_12_9_17';         % same label (variant)
His.pep_seq{no,1} = 'KSHGGKAPR';
His.mod_type{no,1} = '0,pr;1,pr;6,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_12_9_17';         % same label (variant)
His.pep_seq{no,1} = 'ISTGGKAPR';
His.mod_type{no,1} = '0,pr;6,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_13_18_26';        % 18–26 variant set
His.pep_seq{no,1} = 'KQLATKAAR';
His.mod_type{no,1} = '0,pr;1,pr;6,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_13_18_26';        % same label (variant)
His.pep_seq{no,1} = 'KELATKAAR';
His.mod_type{no,1} = '0,pr;1,pr;6,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_13_18_26';        % same label (variant)
His.pep_seq{no,1} = 'TLLATKAAR';
His.mod_type{no,1} = '0,pr;6,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_13_18_26';        % same label (variant)
His.pep_seq{no,1} = 'KQLAPKAAR';
His.mod_type{no,1} = '0,pr;1,pr;6,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_14_27_40';        % 27–40 variant set (key plant biology: H3.1 vs H3.3)
His.pep_seq{no,1} = 'KSAPATGGVKKPHR';         % "AT" motif (H3.1-like)
His.mod_type{no,1} = '0,pr;1,pr;10,pr;11,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_14_27_40';        % same region with the "TT" motif (H3.3-like)
His.pep_seq{no,1} = 'KSAPTTGGVKKPHR';
His.mod_type{no,1} = '0,pr;1,pr;10,pr;11,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_14_27_40';        % same label; Q at N-1 variant
His.pep_seq{no,1} = 'QSAPATGGVKKPHR';
His.mod_type{no,1} = '0,pr;10,pr;11,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_16_53_63';        % duplicates/variants of 53–63 window
His.pep_seq{no,1} = 'KYQKSTELLIR';
His.mod_type{no,1} = '0,pr;1,pr;4,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_16_53_63';        % same label; NR tail variant
His.pep_seq{no,1} = 'KYQKSTELLNR';
His.mod_type{no,1} = '0,pr;1,pr;4,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_17_73_83';        % duplicates/variants of 73–83 window
His.pep_seq{no,1} = 'EIAQDFKTDLR';
His.mod_type{no,1} = '0,pr;7,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_17_73_83';        % same label; Y variant
His.pep_seq{no,1} = 'EIAQDYKTDLR';
His.mod_type{no,1} = '0,pr;7,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_18_117_128';      % duplicates/variants of 117–128 window
His.pep_seq{no,1} = 'VTIMPKDIQLAR';
His.mod_type{no,1} = '0,pr;6,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_18_117_128';      % same label; VQ variant
His.pep_seq{no,1} = 'VTIMPKDVQLAR';
His.mod_type{no,1} = '0,pr;6,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H3_18_117_128';      % same label; EI variant
His.pep_seq{no,1} = 'VTIMPKEIQLAR';
His.mod_type{no,1} = '0,pr;6,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

%------------------H4------------------
no = no + 1;
His.out_filename{no,1} = 'H4_02a_18_23';       % H4 short window variant
His.pep_seq{no,1} = 'HRKVLR';
His.mod_type{no,1} = '0,pr;3,pr;';             % N-term + specific Lys derivatization
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H4_02b_20_35';       % H4 longer window (charge 3)
His.pep_seq{no,1} = 'KVLRDNIQGITKPAIR';
His.mod_type{no,1} = '0,pr;1,pr;12,pr;';
His.pep_ch(no,1) = 3;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H4_02c_20_36';       % H4 longer window (+R)
His.pep_seq{no,1} = 'KVLRDNIQGITKPAIRR';
His.mod_type{no,1} = '0,pr;1,pr;12,pr;';
His.pep_ch(no,1) = 3;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H4_03_24_35';        % H4 internal anchor
His.pep_seq{no,1} = 'DNIQGITKPAIR';
His.mod_type{no,1} = '0,pr;8,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H4_04_40_45';        % H4 short anchor
His.pep_seq{no,1} = 'RGGVKR';
His.mod_type{no,1} = '0,pr;5,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H4_05_68_78';        % H4 internal window
His.pep_seq{no,1} = 'DAVTYTEHAKR';
His.mod_type{no,1} = '0,pr;10,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H4_06_79_92';        % H4 internal window
His.pep_seq{no,1} = 'KTVTAMDVVYALKR';
His.mod_type{no,1} = '0,pr;1,pr;13,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H4_07u_24_102';      % grouped H4 micro-anchors (same label reused)
His.pep_seq{no,1} = 'DNIQGITKPAIRR';
His.mod_type{no,1} = '0,pr;8,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H4_07u_24_102';      % same label (variant)
His.pep_seq{no,1} = 'ISGLIYEETR';
His.mod_type{no,1} = '0,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H4_07u_24_102';      % same label (variant)
His.pep_seq{no,1} = 'GVLKVFLENVIR';
His.mod_type{no,1} = '0,pr;4,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'H4_07u_24_102';      % same label (variant)
His.pep_seq{no,1} = 'TLYGFGG';
His.mod_type{no,1} = '0,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

%------------------H2A------------------
no = no + 1;
His.out_filename{no,1} = 'HH2A_01u_1_7';       % minimal H2A N-term variants (QC anchors; optional)
His.pep_seq{no,1} = 'DNKKSR';
His.mod_type{no,1} = '0,pr;3,pr;4,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'HH2A_01u_1_7';       % same label (variant)
His.pep_seq{no,1} = 'DNKKTR';
His.mod_type{no,1} = '0,pr;3,pr;4,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'HH2A_01u_1_7';       % same label (variant)
His.pep_seq{no,1} = 'DNKKNR';
His.mod_type{no,1} = '0,pr;3,pr;4,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'HH2A_01u_1_7';       % no fixed Lys derivatization here
His.pep_seq{no,1} = 'HLLLAIR';
His.mod_type{no,1} = '0,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'HH2A_01u_1_7';       % variant
His.pep_seq{no,1} = 'HLQLAIR';
His.mod_type{no,1} = '0,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'HH2A_01u_1_7';       % variant
His.pep_seq{no,1} = 'HLCLAIR';
His.mod_type{no,1} = '0,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'HH2A_01u_1_7';       % variant
His.pep_seq{no,1} = 'HIQLAVR';
His.mod_type{no,1} = '0,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

no = no + 1;
His.out_filename{no,1} = 'HH2A_01u_1_7';       % variant
His.pep_seq{no,1} = 'HVLLAVR';
His.mod_type{no,1} = '0,pr;';
His.pep_ch(no,1) = 2;
His.pep_mz(no,1) = calculate_pepmz0(His,no,special);
new_seq = [His.pep_seq{no,1},His.mod_type{no,1}];
His.seq_godel(no,1) = sum((new_seq-'0'+49).*log(2:1+length(new_seq)));

%------------------------------------
function pep_mz = calculate_pepmz0(His,hno,special)
%%
% PURPOSE:
%   Compute theoretical monoisotopic m/z for the peptide indicated by index "hno"
%   using the provided amino-acid mass table selection (via "special") and fixed mods.
%
% INPUTS:
%   His     : struct with peptide fields populated up to row "hno".
%   hno     : integer row index of the peptide to compute.
%   special : struct controlling AA mass table selection (instrument/label mode).
%
% DEPENDENCIES:
%   GetMods()      : returns modification catalogue and masses.
%   get_mod_mass() : returns per-residue delta-mass array for this peptide+mods.
%   Getaamass()    : canonical monoisotopic AA masses.
%   GetaamassH()   : alternative AA masses (only used in specific "special" mode).

Mods = GetMods();
if 4==special.nsource && (1==special.nsubtype || 3==special.nsubtype)
    % Use alternative AA mass table (e.g., special instrument/label mode)
    aamass = GetaamassH();
else
    % Default AA mass table
    aamass = Getaamass();
end;
element = [12 1.0078246 14.0030732 15.9949141 31.972070];% elemental masses: C, H, N, O, S
mH2O = element(2)*2 + element(4);                         % add H2O for peptide mass (closed ends)
pmass = 1.007276;                                         % proton mass for charge handling

c_seq = His.pep_seq{hno};                                 % current peptide sequence
idx = c_seq-'A'+1;                                        % map letters to rows in AA mass table
residuemass = aamass(idx,1)';                             % row vector of residue masses
c_mod = His.mod_type{hno};                                % fixed modifications "pos,type;" grammar
deltam = get_mod_mass(c_seq,c_mod,Mods);                  % per-position delta masses from mods
% peptide+modification
residuemass_new = residuemass + deltam;                   % apply fixed mods
Mr = sum(residuemass_new)+mH2O;                           % neutral peptide mass
c_ch = His.pep_ch(hno);                                   % target charge
pep_mz = (Mr+c_ch*pmass)/c_ch;                            % convert to m/z at charge c_ch

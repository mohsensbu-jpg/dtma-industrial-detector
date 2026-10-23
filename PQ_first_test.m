% =========================================================================
% PQ_first_test.m
% Power-quality impulsive-transient detection with Mohsen Karimi's
% D-TMA band-pass + local decision rule (ideas: Mohsen Karimi).
% Step 1 of transferring the ECG/QRS detector to an industrial domain.
%
% - Downloads ONE real voltage recording (Mendeley Data 5jp4r7nnpw v2,
%   "BidyutPi" 10 kHz waveforms, 330 W resistive load, ~30 MB) into
%   MATLAB Drive (the download happens on the MATLAB Online server).
% - Exp A : real grid voltage as background + injected transients
%           + REAL recorded noise (residual of the same recording) at
%           controlled SNR (24/18/12/6/0 dB) -> "decisive band" protocol.
%           A1 = noise as recorded, A2 = same noise with slow level changes.
% - Exp B : fully synthetic IEEE-1159-style grid, fs = 50 kHz,
%           non-stationary noise (runs even if download fails).
% - Part C: detector on the raw real recording (no injection)
%           -> real false-alarm rate on a clean grid.
% Text output only. No toolbox functions.
% Parameters are set physically (W, 10-cycle window, minimum interval)
% and are NOT tuned on these data.
% =========================================================================
clear; clc; tic;
fprintf('==== PQ_first_test  (%s) ====\n', datestr(now));
SNRs = [24 18 12 6 0];
nRec = 20;            % records per SNR level
A0   = 0.2;           % transient amplitude [p.u. of fundamental peak]

%% ---------------- Part 0 : download + read real data -------------------
fn  = 'V_I_Waveforms_330W_Resistive_Load.csv';
url = ['https://data.mendeley.com/public-files/datasets/5jp4r7nnpw/files/' ...
       'e37d5ad6-dd42-48ce-99e6-8ef6b5a6e3d8/file_downloaded'];
haveReal = false;
if ~isfile(fn)
    try
        fprintf('Downloading %s (~30 MB) ...\n', fn);
        websave(fn, url, weboptions('Timeout', 600));
        fprintf('Download OK.\n');
    catch ME
        fprintf('DOWNLOAD FAILED: %s\n', ME.message);
    end
end
if isfile(fn)
    try
        T = readtable(fn, 'VariableNamingRule', 'preserve');
        names = T.Properties.VariableNames;
        fprintf('\nFile: %s   rows = %d\nColumns: %s\nFirst rows:\n', ...
            fn, height(T), strjoin(names, ' | '));
        disp(T(1:3,:));
        isnum = false(1, numel(names));
        for c = 1:numel(names), isnum(c) = isnumeric(T.(names{c})); end
        ln = lower(names);
        vi = find(isnum & contains(ln, 'volt'), 1);
        if isempty(vi), vi = find(isnum & startsWith(ln, 'v'), 1); end
        ti = find(isnum & (contains(ln, 'time') | strcmp(ln, 't')), 1);
        if isempty(vi)                      % fallback: largest-RMS numeric column
            best = -1;
            for c = find(isnum)
                if c == ti, continue; end
                r = std(double(T.(names{c})), 'omitnan');
                if r > best, best = r; vi = c; end
            end
        end
        v = double(T.(names{vi})); v = v(:); v(isnan(v)) = [];
        v = v - mean(v);
        Nv = numel(v);
        [tc, nUp] = upcross(v);                 % up-crossing positions (fractional samples)
        cps = (nUp - 1) / (tc(end) - tc(1));    % cycles per sample
        fs = NaN;
        if ~isempty(ti)
            dt = median(diff(double(T.(names{ti}))));
            for sc = [1 1e-3 1e-6]
                fsc = 1 / (dt * sc);
                if cps * fsc > 40 && cps * fsc < 70, fs = fsc; break; end
            end
        end
        if isnan(fs)
            fs = 50 / cps;
            fprintf('NOTE: sampling rate not found in time column; set from f0 = 50 Hz.\n');
        end
        f0 = cps * fs;
        A1 = sqrt(2) * sqrt(mean(v.^2));
        vpu = v / A1;                            % per-unit of fundamental peak (approx.)
        if isempty(ti), tname = '-'; else, tname = names{ti}; end
        fprintf('Voltage column: "%s"   time column: "%s"\n', names{vi}, tname);
        fprintf('fs = %.1f Hz   N = %d (%.1f s)   f0 = %.3f Hz   Vpeak = %.2f (raw units)\n', ...
            fs, Nv, Nv/fs, f0, A1);
        % ---- real noise = residual after block-wise harmonic LS fit ----
        nzr = harmonicResidual(vpu, fs, tc);
        sig_n = sqrt(mean(nzr.^2));
        fprintf('Real residual noise (aperiodic part): RMS = %.5f p.u.  -> intrinsic SNR for A0=%.2f: %.1f dB\n', ...
            sig_n, A0, 20*log10(A0/sig_n));
        haveReal = true;
    catch ME
        fprintf('READING/PARSING FAILED: %s\n', ME.message);
    end
end

%% ---------------- Exp A : real background + real noise -----------------
if haveReal
    Wsec = 0.5e-3;  tau1 = 100e-6;  tau2 = 400e-6;   % ms-scale impulsive transient
    [sh, pk] = shape(fs, tau1, tau2);
    Lseg = round(4*fs); if 3*Lseg > numel(nzr), Lseg = floor(numel(nzr)/3); end
    for modn = [false true]
        rng(11);
        recFun = @(snr) recReal(vpu, nzr, fs, Lseg, snr, A0, sh, pk, modn);
        tag = sprintf('EXP A%d: REAL grid voltage + REAL noise (%s), fs=%.0f Hz', ...
            1+modn, ternary(modn, 'with slow level changes 0.4x..2x', 'as recorded'), fs);
        runExp(tag, fs, f0, Wsec, nRec, SNRs, recFun, tau1, tau2, A0);
    end
else
    fprintf('\n*** Exp A skipped (no real data). ***\n');
end

%% ---------------- Exp B : synthetic grid, 50 kHz ------------------------
fsB = 50000; f0B = 50; WsecB = 0.2e-3; t1B = 20e-6; t2B = 150e-6;
[shB, pkB] = shape(fsB, t1B, t2B);
rng(22);
recFunB = @(snr) recSynth(fsB, round(4*fsB), snr, A0, shB, pkB);
runExp('EXP B: SYNTHETIC grid (harmonics 3,5,7,11, flicker, f-drift) + non-stationary noise, fs=50 kHz', ...
    fsB, f0B, WsecB, nRec, SNRs, recFunB, t1B, t2B, A0);

%% ---------------- Part C : raw real recording, no injection -------------
if haveReal
    P = params(fs, 0.5e-3, f0);
    [d1, sc1] = detMain(vpu, P);
    d2 = detECG(vpu, P);
    fprintf('\n==== PART C: raw real recording (%.1f s), nothing injected ====\n', Nv/fs);
    fprintf('Main detector : %d detections (%.2f /s)\n', numel(d1), numel(d1)/(Nv/fs));
    fprintf('ECG rule      : %d detections (%.2f /s)\n', numel(d2), numel(d2)/(Nv/fs));
    if ~isempty(d1)
        [~, o] = sort(sc1(d1), 'descend'); o = o(1:min(8, numel(o)));
        fprintf('Strongest main detections (time s, e/s):\n');
        for j = o(:)', fprintf('   t = %9.4f   e/s = %6.2f\n', d1(j)/fs, sc1(d1(j))); end
    end
end
fprintf('\nTotal run time: %.1f s\n==== END ====\n', toc);

% =========================================================================
%                              FUNCTIONS
% =========================================================================
function P = params(fs, Wsec, f0)
P.fs = fs; P.f0 = f0;
P.W  = max(3, round(Wsec*fs));
P.w1 = max(3, round(3.3*P.W));
P.w2 = max(3, round(0.06*P.W));
P.L  = round(10/f0*fs);          % 10 fundamental cycles (IEC 61000-4-30 window)
P.Tmin = round(2e-3*fs);         % physical minimum interval between transients
P.tol  = max(2, round(0.5e-3*fs));
P.N0 = round(fs/f0);             % one cycle (for cycle-difference comparator)
end

function runExp(tag, fs, f0, Wsec, nRec, SNRs, recFun, tau1, tau2, A0)
P = params(fs, Wsec, f0);
M = {'Karimi D-TMA + local robust rule', 'Karimi ECG rule unchanged', ...
     'D-TMA + GLOBAL thr (ablation)', 'db4 wavelet d1, universal thr', ...
     'Cycle difference, universal thr'};
nM = numel(M);
fprintf('\n=====================================================================\n');
fprintf('%s\n', tag);
fprintf('transient: double-exp tau1=%.0f us tau2=%.0f us, amp %.2f p.u. x U(0.6,1.4), random sign\n', ...
    tau1*1e6, tau2*1e6, A0);
fprintf('W=%d smp (%.2f ms)  w1=%d  w2=%d  L=%d (10 cycles)  Tmin=%d  tol=+-%d smp (%.2f ms)\n', ...
    P.W, P.W/fs*1e3, P.w1, P.w2, P.L, P.Tmin, P.tol, P.tol/fs*1e3);
fprintf('local rule: thr = sqrt(2 ln L) * 1.4826*movMAD = %.2f * scale;  %d records x 4 s per SNR\n', ...
    sqrt(2*log(P.L)), nRec);
fprintf('SNR = 20 log10(A0 / RMS(added noise))\n');
fprintf('=====================================================================\n');
winSig = cell(1, nM);
for si = 1:numel(SNRs)
    snr = SNRs(si);
    C = zeros(nRec, nM, 3); errs = []; dur = 0;
    for r = 1:nRec
        [x, ev] = recFun(snr);
        dur = dur + numel(x)/fs;
        D = cell(1, nM);
        D{1} = detMain(x, P);  D{2} = detECG(x, P);  D{3} = detGlobal(x, P);
        D{4} = detWavelet(x, P);  D{5} = detCycle(x, P);
        for m = 1:nM
            [tp, fp, fn, e] = matchEv(D{m}, ev, P.tol);
            C(r, m, :) = [tp fp fn];
            if m == 1, errs = [errs; e(:)]; end %#ok<AGROW>
        end
    end
    fprintf('\n--- SNR = %2d dB  (%d true events) ---\n', snr, sum(C(:,1,1)) + sum(C(:,1,3)));
    for m = 1:nM
        tp = sum(C(:,m,1)); fp = sum(C(:,m,2)); fn = sum(C(:,m,3));
        F1 = 2*tp/max(1, 2*tp+fp+fn); Se = tp/max(1, tp+fn); PPV = tp/max(1, tp+fp);
        fprintf('%-34s F1=%.3f Se=%.3f PPV=%.3f FP/s=%6.2f', M{m}, F1, Se, PPV, fp/dur);
        if m == 1
            if ~isempty(errs)
                off = myq(errs, 0.5);
                fprintf('  offset=%+.0f us  |err| q50/q95=%.0f/%.0f us\n', off/fs*1e6, ...
                    myq(abs(errs-off), 0.5)/fs*1e6, myq(abs(errs-off), 0.95)/fs*1e6);
            else
                fprintf('\n');
            end
        else
            [dm, lo, hi] = bootDiff(C(:,1,:), C(:,m,:), 2000);
            star = '';
            if lo > 0, star = '  * main better'; winSig{m}(end+1) = snr; end
            if hi < 0, star = '  * this better'; end
            fprintf('  dF1(main-this)=%+.3f [%+.3f %+.3f]%s\n', dm, lo, hi, star);
        end
    end
end
fprintf('\nSummary (main significantly better, 95%% bootstrap CI > 0):\n');
for m = 2:nM
    fprintf('  vs %-33s at SNR: %s\n', M{m}, mat2str(winSig{m}));
end
end

% ---------------- detectors ----------------------------------------------
function [idx, score] = detMain(x, P)
xp = dtma(x, P);
mu = movmedian(xp, P.L);
e  = abs(xp - mu);
s  = max(1.4826 * movmedian(e, P.L), eps);
k  = sqrt(2*log(P.L));
score = e ./ s;
idx = lmax(e, k*s);
idx = mergeg(idx, score, P.w1);
idx = mergeg(idx, score, P.Tmin);
end

function idx = detECG(x, P)      % his QRS rule, only window/merge rescaled
xp = dtma(x, P);
d  = xp - movmedian(xp, P.L);
s  = max(movstd(d, P.L), eps);
e  = abs(d); score = e ./ s;
idx = lmax(e, 1.75*s);
idx = mergeg(idx, score, P.w1);
idx = mergeg(idx, score, P.Tmin);
end

function idx = detGlobal(x, P)   % ablation: same filter, global robust scale
xp = dtma(x, P);
e  = abs(xp - median(xp));
s  = max(1.4826 * median(e), eps);
idx = lmax(e, sqrt(2*log(P.L))*s);
idx = mergeg(idx, e, P.w1);
idx = mergeg(idx, e, P.Tmin);
end

function idx = detWavelet(x, P)  % classic PQ: db4 level-1 detail (undecimated)
g = [-0.0105974018 -0.0328830117 0.0308413818 0.1870348117 ...
     -0.0279837694 -0.6308807679 0.7148465706 -0.2303778133];
d = conv(x(:), g(:), 'same');
e = abs(d);
thr = median(e)/0.6745 * sqrt(2*log(numel(x)));
idx = lmax(e, thr);
idx = mergeg(idx, e, P.Tmin);
end

function idx = detCycle(x, P)    % classic PQ: point-by-point cycle comparison
x = x(:); r = zeros(size(x));
r(P.N0+1:end) = x(P.N0+1:end) - x(1:end-P.N0);
e = abs(r);
thr = median(e(P.N0+1:end))/0.6745 * sqrt(2*log(numel(x)));
idx = lmax(e, thr);
idx = mergeg(idx, e, P.Tmin);
end

% ---------------- filter --------------------------------------------------
function xp = dtma(x, P)
x  = x(:);
xp = tma(x - tma(x, P.w1), P.w2);
end

function y = tma(x, w)          % two cascaded length-K boxcars, zero phase
K = ceil(w/2);
y = movmean(x, K);
if mod(K, 2) == 0
    y = movmean(y, [K/2-1 K/2]);  % cancels the half-sample shift of the first
else
    y = movmean(y, K);
end
end

% ---------------- helpers -------------------------------------------------
function idx = lmax(e, thr)
e = e(:); n = numel(e);
if isscalar(thr), thr = thr*ones(n, 1); else, thr = thr(:); end
i = (2:n-1)';
m = e(i) > e(i-1) & e(i) >= e(i+1) & e(i) > thr(i);
idx = i(m);
end

function out = mergeg(idx, score, gap)
out = zeros(numel(idx), 1); c = 0;
for j = 1:numel(idx)
    i = idx(j);
    if c > 0 && i - out(c) < gap
        if score(i) > score(out(c)), out(c) = i; end
    else
        c = c + 1; out(c) = i;
    end
end
out = out(1:c);
end

function [tp, fp, fn, err] = matchEv(det, ev, tol)
det = det(:); used = false(size(det)); tp = 0; err = [];
for j = 1:numel(ev)
    if isempty(det), break; end
    d = abs(det - ev(j)); d(used) = inf;
    [dm, k] = min(d);
    if dm <= tol
        used(k) = true; tp = tp + 1; err(end+1, 1) = det(k) - ev(j); %#ok<AGROW>
    end
end
fp = numel(det) - tp; fn = numel(ev) - tp;
end

function [dm, lo, hi] = bootDiff(Ca, Cb, B)
Ca = reshape(Ca, [], 3); Cb = reshape(Cb, [], 3); n = size(Ca, 1);
f = @(C) 2*sum(C(:,1)) / max(1, 2*sum(C(:,1)) + sum(C(:,2)) + sum(C(:,3)));
dm = f(Ca) - f(Cb); d = zeros(B, 1);
for b = 1:B
    ii = randi(n, n, 1);
    d(b) = f(Ca(ii,:)) - f(Cb(ii,:));
end
lo = myq(d, 0.025); hi = myq(d, 0.975);
end

function q = myq(x, p)          % own quantile (no toolbox)
x = sort(x(:)); n = numel(x);
if n == 0, q = NaN; return; end
pos = 1 + (n-1)*p; lo = floor(pos); hi = ceil(pos);
q = x(lo) + (pos-lo)*(x(hi) - x(lo));
end

function [sh, pk] = shape(fs, tau1, tau2)
tt = (0:round(10*tau2*fs))'/fs;
sh = exp(-tt/tau2) - exp(-tt/tau1);
sh = sh / max(sh);
[~, pk] = max(sh);
end

function [x, ev] = inject(bg, nz, fs, snr, A0, sh, pk, modn)
N = numel(bg); ev = []; p = round(0.02*fs);
while true
    p = p + round((0.005 + 0.115*rand)*fs);
    if p > N - round(0.02*fs) - numel(sh), break; end
    ev(end+1, 1) = p; %#ok<AGROW>
end
s = zeros(N, 1);
for j = 1:numel(ev)
    a  = (2*(rand > 0.5) - 1) * (0.6 + 0.8*rand) * A0;
    ii = (ev(j)-pk+1):(ev(j)-pk+numel(sh));
    s(ii) = s(ii) + a*sh;
end
n = nz(:);
if modn
    g = interp1(linspace(1, N, 9), 0.4 + 1.6*rand(1, 9), (1:N)');
    n = n .* g;
end
n = n / sqrt(mean(n.^2)) * A0 / 10^(snr/20);
x = bg(:) + s + n;
end

function [x, ev] = recReal(vpu, nzr, fs, Lseg, snr, A0, sh, pk, modn)
Nv = numel(nzr);                 % nzr is slightly shorter than vpu
s1 = randi(Nv - Lseg + 1);
s2 = s1;
while abs(s2 - s1) < Lseg, s2 = randi(Nv - Lseg + 1); end
[x, ev] = inject(vpu(s1:s1+Lseg-1), nzr(s2:s2+Lseg-1), fs, snr, A0, sh, pk, modn);
end

function [x, ev] = recSynth(fs, N, snr, A0, sh, pk)
t   = (0:N-1)'/fs;
f   = 50 + 0.1*sin(2*pi*0.3*t + 2*pi*rand);
ph  = 2*pi*cumsum(f)/fs;
amp = 1 + 0.02*sin(2*pi*8*t);
bg  = amp.*sin(ph) + 0.03*sin(3*ph+1) + 0.04*sin(5*ph+2) + 0.025*sin(7*ph) + 0.01*sin(11*ph);
nz  = randn(N, 1) + movmean(randn(N, 1), 5)*sqrt(5);
[x, ev] = inject(bg, nz, fs, snr, A0, sh, pk, true);
end

function [tc, nUp] = upcross(x)  % up-crossings with hysteresis, linear interpolation
a = 0.1*max(abs(x)); st = 0; tc = zeros(ceil(numel(x)/2), 1); nUp = 0; last = 1;
for i = 2:numel(x)
    if x(i-1) <= 0 && x(i) > 0, last = i; end
    if st <= 0 && x(i) > a
        if st < 0
            k = last;
            nUp = nUp + 1;
            tc(nUp) = (k-1) + (-x(k-1))/(x(k) - x(k-1));
        end
        st = 1;
    elseif st >= 0 && x(i) < -a
        st = -1;
    end
end
tc = tc(1:nUp);
end

function nz = harmonicResidual(v, fs, tc)
% block-wise (0.2 s) LS fit of DC + all harmonics below 0.45 fs, with a
% block fundamental from the zero crossings -> residual = aperiodic part
N = numel(v); B = round(0.2*fs); nz = zeros(N, 1);
f0g = (numel(tc)-1)/(tc(end)-tc(1))*fs;
for b = 1:floor(N/B)
    ii = ((b-1)*B + 1 : b*B)';
    c  = tc(tc >= ii(1)-1 & tc <= ii(end)-1);
    if numel(c) >= 3, f0 = (numel(c)-1)/(c(end)-c(1))*fs; else, f0 = f0g; end
    t  = (0:B-1)'/fs;
    nh = min(floor(0.45*fs/f0), 100);
    H  = ones(B, 1 + 2*nh);
    for h = 1:nh
        H(:, 2*h)   = cos(2*pi*h*f0*t);
        H(:, 2*h+1) = sin(2*pi*h*f0*t);
    end
    cf = H \ v(ii);
    nz(ii) = v(ii) - H*cf;
end
nb = floor(N/B)*B;
nz = nz(1:nb);                   % drop incomplete last block
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end

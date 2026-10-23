% =========================================================================
% PQ_four_loads.m -- real-data test of Mohsen Karimi's D-TMA + local rule
% Data: Mendeley Data 5jp4r7nnpw v2 (BidyutPi, Raspberry-Pi meter, 10 kHz),
%       4 recordings: resistive, induction motor, ceiling fan, SMPS charger.
% Design choices (lessons from the first test):
%   * no "clean segment" search any more: recorder gaps (+-20 ms) are
%     MASKED - no events are injected there and detections there are ignored
%   * SHORT CALIBRATION (his idea: online use after a short calibration):
%     first 15 s of each recording = calibration only. The threshold factor k
%     is the smallest k giving <= 0.1 false alarms/s there (nothing injected).
%     The same calibration is given to the db4 wavelet (fair comparison).
%     All tests use only the remaining 45 s.
%   * reference for the +/- marks = main_cal (calibrated main method)
% Text output only, no toolbox. Downloads files 3-4 (~60 MB) on the server.
% =========================================================================
clear; clc; tic;
fprintf('==== PQ_four_loads  (%s) ====\n', datestr(now));
base  = 'https://data.mendeley.com/public-files/datasets/5jp4r7nnpw/files/';
files = {'V_I_Waveforms_330W_Resistive_Load.csv',          'e37d5ad6-dd42-48ce-99e6-8ef6b5a6e3d8', 'resistive 330 W';
         'V_I_Waveforms_0.5HP_Induction_Motor_No_Load.csv', 'c3fdd71a-3b72-400b-92e7-3592b04c414a', 'induction motor';
         'V_I_Waveforms_50W_Ceiling_Fan_Aged.csv',          'ecc9b904-7879-4892-84d0-d78fc7fef9a6', 'ceiling fan 50 W';
         'V_I_Waveforms_SMPS_Mobile_Charger.csv',           '795ae74d-afde-431d-abe6-3c50b0002096', 'SMPS charger'};
SNRs = [24 18 12 6 0];  nRec = 20;  A0 = 0.2;  Tcal = 15;  FAtarget = 0.1;
Wsec = 0.5e-3;  tau1 = 100e-6;  tau2 = 400e-6;
MN = {'main', 'ECG', 'global', 'db4', 'cycle', 'mainCAL', 'db4CAL'};  nM = 7;  REF = 6;
nF = size(files,1);
FA = nan(nF, nM);  F12 = nan(nF, nM);  F18 = nan(nF, nM);  SNRi = nan(nF,1);  KC = nan(nF,2);

for f = 1:nF
    fn = files{f,1};
    fprintf('\n#####################################################################\n');
    fprintf('FILE %d: %s  (%s)\n', f, files{f,3}, fn);
    if ~isfile(fn)
        try
            fprintf('downloading (~30 MB) ...\n');
            websave(fn, [base files{f,2} '/file_downloaded'], weboptions('Timeout', 600));
        catch ME
            fprintf('DOWNLOAD FAILED: %s\n', ME.message); continue;
        end
    end
    try
        [x, fs, f0, bad, ng] = loadUniform(fn);
    catch ME
        fprintf('READ FAILED: %s\n', ME.message); continue;
    end
    N = numel(x);  P = params(fs, Wsec, f0);
    [tc, ~] = upcross(x);
    nz = harmonicResidual(x, fs, tc);
    bad(numel(nz)+1:end) = true;  nz(end+1:N) = 0;
    sig_i = sqrt(mean(nz(~bad).^2));  SNRi(f) = 20*log10(A0/sig_i);
    fprintf('fs=%d Hz  N=%d (%.1f s)  f0=%.3f Hz  recorder gaps=%d  masked=%.2f %%\n', ...
        fs, N, N/fs, f0, ng, 100*mean(bad));
    fprintf('intrinsic noise RMS = %.4f p.u. -> intrinsic SNR (A0=%.2f) = %.1f dB\n', sig_i, A0, SNRi(f));

    % ---------- calibration on the first Tcal seconds ----------
    ic = 1:round(Tcal*fs);  it = (round(Tcal*fs)+1):N;
    [kM, kW, faM, faW] = calibrate(x(ic), bad(ic), P, FAtarget);
    KC(f,:) = [kM kW];
    Pc = P;  Pc.k = kM;  Pc.kw = kW;
    fprintf('calibration (first %d s, target <= %.2f FA/s): main k = %.2f (default %.2f, FA/s %.3f) | db4 k = %.2f (default %.2f, FA/s %.3f)\n', ...
        Tcal, FAtarget, kM, P.k, faM, kW, sqrt(2*log(4*fs)), faW);

    % ---------- Part 1: false alarms on the TEST part, nothing injected ----------
    xt = x(it);  bt = bad(it);  Tclean = sum(~bt)/fs;
    D = runAll(xt, P, Pc);
    fprintf('Part 1  false alarms/s on real voltage (test part, %.1f s clean):\n  ', Tclean);
    for m = 1:nM
        d = D{m}; d = d(~bt(d)); FA(f,m) = numel(d)/Tclean;
        fprintf('%s=%.2f  ', MN{m}, FA(f,m));
    end
    fprintf('\n');

    % ---------- Part 2: injection, total-noise SNR ----------
    [sh, pk] = shape(fs, tau1, tau2);
    Lseg = round(4*fs);  pool = nz(it);  pool = pool(~bt);
    rng(100 + f);
    fprintf('Part 2  F1; marks vs mainCAL: + mainCAL significantly better, - worse, = n.s.\n');
    for si = 1:numel(SNRs)
        sig_t = A0/10^(SNRs(si)/20);
        if sig_t <= sig_i
            fprintf('  SNR %2d: not reachable (intrinsic noise already higher)\n', SNRs(si)); continue;
        end
        sig_a = sqrt(sig_t^2 - sig_i^2);
        C = zeros(nRec, nM, 3);
        for r = 1:nRec
            s1 = randi(numel(xt) - Lseg + 1);  ii = s1:s1+Lseg-1;
            s2 = randi(numel(pool) - Lseg + 1);
            [xr, ev] = inject3(xt(ii), bt(ii), pool(s2:s2+Lseg-1), fs, sig_a, A0, sh, pk);
            Dr = runAll(xr, P, Pc);
            for m = 1:nM
                d = Dr{m}; d = d(~bt(ii(d)));
                [tp, fp, fnn] = matchEv(d, ev, P.tol);
                C(r, m, :) = [tp fp fnn];
            end
        end
        fprintf('  SNR %2d |', SNRs(si));
        for m = 1:nM
            tp = sum(C(:,m,1)); F1 = 2*tp/max(1, 2*tp+sum(C(:,m,2))+sum(C(:,m,3)));
            if SNRs(si) == 12, F12(f,m) = F1; end
            if SNRs(si) == 18, F18(f,m) = F1; end
            if m == REF
                fprintf(' %s %.3f (Se %.2f PPV %.2f)', MN{m}, F1, tp/max(1,tp+sum(C(:,m,3))), tp/max(1,tp+sum(C(:,m,2))));
            else
                [~, lo, hi] = bootDiff(C(:,REF,:), C(:,m,:), 2000);
                mk = '='; if lo > 0, mk = '+'; elseif hi < 0, mk = '-'; end
                fprintf(' %s %.3f%s', MN{m}, F1, mk);
            end
        end
        fprintf('\n');
    end
end

fprintf('\n#####################################################################\nSUMMARY\n');
fprintf('calibrated k:  '); for f = 1:nF, fprintf('%s: main %.2f, db4 %.2f;  ', files{f,3}, KC(f,1), KC(f,2)); end
fprintf('\n\n%-17s %5s | false alarms/s on real voltage (test part)\n', 'file', 'SNRi');
fprintf('%-17s %5s |', '', ''); fprintf(' %7s', MN{:}); fprintf('\n');
for f = 1:nF, fprintf('%-17s %5.1f |', files{f,3}, SNRi(f)); fprintf(' %7.2f', FA(f,:)); fprintf('\n'); end
fprintf('\n%-17s %5s | F1 at total SNR 18 dB\n', 'file', '');
for f = 1:nF, fprintf('%-17s %5s |', files{f,3}, ''); fprintf(' %7.3f', F18(f,:)); fprintf('\n'); end
fprintf('\n%-17s %5s | F1 at total SNR 12 dB\n', 'file', '');
for f = 1:nF, fprintf('%-17s %5s |', files{f,3}, ''); fprintf(' %7.3f', F12(f,:)); fprintf('\n'); end
fprintf('\nTotal run time: %.1f s\n==== END ====\n', toc);

% =========================================================================
%                              FUNCTIONS
% =========================================================================
function D = runAll(x, P, Pc)
D = {detMain(x, P), detECG(x, P), detGlobal(x, P), detWavelet(x, P), detCycle(x, P), ...
     detMain(x, Pc), detWavelet(x, Pc)};
end

function [kM, kW, faM, faW] = calibrate(x, bad, P, target)
T = sum(~bad)/P.fs;  ks = 3:0.25:20;
kM = ks(end); faM = NaN; kW = ks(end); faW = NaN;
for kk = ks
    Q = P; Q.k = kk; d = detMain(x, Q); fa = sum(~bad(d))/T;
    if fa <= target, kM = kk; faM = fa; break; end
end
for kk = ks
    Q = P; Q.kw = kk; d = detWavelet(x, Q); fa = sum(~bad(d))/T;
    if fa <= target, kW = kk; faW = fa; break; end
end
end

function [x, ev] = inject3(bg, bseg, nz, fs, sig_a, A0, sh, pk)
N = numel(bg); ev = []; p = round(0.02*fs); mg = round(0.01*fs);
while true
    p = p + round((0.005 + 0.115*rand)*fs);
    if p > N - round(0.02*fs) - numel(sh), break; end
    if ~any(bseg(max(1,p-mg) : min(N, p+mg+numel(sh)))), ev(end+1, 1) = p; end %#ok<AGROW>
end
s = zeros(N, 1);
for j = 1:numel(ev)
    a  = (2*(rand > 0.5) - 1) * (0.6 + 0.8*rand) * A0;
    ii = (ev(j)-pk+1):(ev(j)-pk+numel(sh));
    s(ii) = s(ii) + a*sh;
end
n = nz(:); n = n/max(sqrt(mean(n.^2)), eps)*sig_a;
x = bg(:) + s + n;
end
function [x, fs, f0, bad, ng] = loadUniform(fn)
T = readtable(fn, 'VariableNamingRule', 'preserve');
t = double(T.time_s(:)); v = double(T.voltage_V(:));
ok = ~isnan(t) & ~isnan(v); t = t(ok); v = v(ok);
[t, ia] = unique(t); v = v(ia);
dt = diff(t); md = median(dt); fs = round(1/md);
tu = (t(1) : 1/fs : t(end))';
x  = interp1(t, v, tu, 'linear');
g  = find(dt > 2*md); ng = numel(g);
bad = false(numel(tu), 1); pad = 0.02;
for k = 1:ng
    i1 = max(1, floor((t(g(k)) - pad - t(1))*fs) + 1);
    i2 = min(numel(tu), ceil((t(g(k)+1) + pad - t(1))*fs) + 1);
    bad(i1:i2) = true;
end
x = x - mean(x);
[tc, nUp] = upcross(x); f0 = (nUp-1)/(tc(end)-tc(1))*fs;
x = x/(sqrt(2)*sqrt(mean(x.^2)));
end

function P = params(fs, Wsec, f0)
P.fs = fs; P.f0 = f0;
P.W  = max(3, round(Wsec*fs));
P.w1 = max(3, round(3.3*P.W));
P.w2 = max(3, round(0.06*P.W));
P.L  = round(10/f0*fs);          % 10 fundamental cycles (IEC 61000-4-30 window)
P.Tmin = round(2e-3*fs);         % physical minimum interval between transients
P.tol  = max(2, round(0.5e-3*fs));
P.N0 = round(fs/f0);             % one cycle (for cycle-difference comparator)
P.k  = sqrt(2*log(P.L));        % threshold factor (default, not calibrated)
P.kw = NaN;                      % NaN -> universal threshold for db4
end

function [idx, score] = detMain(x, P)
xp = dtma(x, P);
mu = movmedian(xp, P.L);
e  = abs(xp - mu);
s  = max(1.4826 * movmedian(e, P.L), eps);
k  = P.k;
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
if isnan(P.kw), kw = sqrt(2*log(numel(x))); else, kw = P.kw; end
thr = median(e)/0.6745 * kw;
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

function [tc, nUp] = upcross(x)  % up-crossings with hysteresis, linear interpolation
a = 0.1*max(abs(x)); st = 0; tc = zeros(ceil(numel(x)/2), 1); nUp = 0; last = 2;
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


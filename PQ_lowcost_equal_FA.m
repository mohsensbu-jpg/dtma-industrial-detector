% =========================================================================
% PQ_lowcost_equal_FA.m
% Portfolio test: does Mohsen Karimi's D-TMA + local rule work on a
% LOW-COST meter (Raspberry-Pi, 10 kHz, Mendeley 5jp4r7nnpw v2) AS IT IS?
% - NO noise is added: only the meter's own real noise.
% - Transients of several amplitudes (0.5 ... 0.1 p.u.) are injected into
%   the real voltage (test part = last 45 s of each recording).
% - Short calibration (his idea): first 15 s, nothing injected. For each
%   false-alarm budget (1, 0.3, 0.1, 0.03 per s) the smallest threshold
%   factor k meeting it is chosen - for D-TMA AND for the db4 wavelet.
% - Result: sensitivity (Se) of both methods at EQUAL false-alarm rate,
%   + measured false alarms on the test part, + block-bootstrap CI (5 s blocks).
% Text output only, no toolbox. Uses the 4 CSV files already in MATLAB Drive.
% =========================================================================
clear; clc; tic;
fprintf('==== PQ_lowcost_equal_FA  (%s) ====\n', datestr(now));
files = {'V_I_Waveforms_330W_Resistive_Load.csv',          'resistive 330 W';
         'V_I_Waveforms_0.5HP_Induction_Motor_No_Load.csv', 'induction motor';
         'V_I_Waveforms_50W_Ceiling_Fan_Aged.csv',          'ceiling fan 50 W';
         'V_I_Waveforms_SMPS_Mobile_Charger.csv',           'SMPS charger'};
Amps = [0.5 0.3 0.2 0.15 0.1];  FAt = [1 0.3 0.1 0.03];  Tcal = 15;  nRep = 3;  Tblk = 5;
Wsec = 0.5e-3;  tau1 = 100e-6;  tau2 = 400e-6;
nA = numel(Amps); nT = numel(FAt); nF = size(files,1);
% pooled block counts: {a,t} -> [tpMain fnMain tpWav fnWav] per block
POOL = cell(nA, nT);
SNRtab = nan(nF, nA);

for f = 1:nF
    fn = files{f,1};
    fprintf('\n#####################################################################\n');
    fprintf('FILE %d: %s\n', f, files{f,2});
    if ~isfile(fn), fprintf('file missing - run PQ_four_loads.m first\n'); continue; end
    [x, fs, f0, bad, ng] = loadUniform(fn);
    N = numel(x);  P = params(fs, Wsec, f0);
    [tc, ~] = upcross(x);  nz = harmonicResidual(x, fs, tc);
    bad(numel(nz)+1:end) = true;  nz(end+1:N) = 0;
    sig_i = sqrt(mean(nz(~bad).^2));
    SNRtab(f,:) = 20*log10(Amps/sig_i);
    fprintf('f0=%.3f Hz  gaps=%d  masked=%.1f %%  noise RMS=%.4f p.u.\n', f0, ng, 100*mean(bad), sig_i);

    ic = 1:round(Tcal*fs);  it = (round(Tcal*fs)+1):N;
    % ---- calibration: one score computation, many thresholds ----
    [eM, sM] = scoreMain(x(ic), P);  [eW, sW] = scoreWav(x(ic));
    Tc = sum(~bad(ic))/fs;  ks = 3:0.05:25;
    kM = nan(1,nT); kW = nan(1,nT);
    for kk = ks                       % thresholds rise -> false alarms fall
        if any(isnan(kM))
            fa = sum(~bad(ic(pickMain(eM, sM, kk, P))))/Tc;
            kM(isnan(kM) & fa <= FAt) = kk;
        end
        if any(isnan(kW))
            fa = sum(~bad(ic(pickWav(eW, sW, kk, P))))/Tc;
            kW(isnan(kW) & fa <= FAt) = kk;
        end
        if ~any(isnan([kM kW])), break; end
    end
    kM(isnan(kM)) = ks(end); kW(isnan(kW)) = ks(end);
    % ---- false alarms on the test part, nothing injected ----
    xt = x(it); bt = bad(it); Tt = sum(~bt)/fs;
    [eM, sM] = scoreMain(xt, P);  [eW, sW] = scoreWav(xt);
    fprintf('budget FA/s :'); fprintf('   %6.2f       ', FAt); fprintf('\n');
    fprintf('k   main/db4:'); fprintf('   %5.2f/%-5.2f  ', [kM; kW]); fprintf('\n');
    fprintf('test FA/s   :');
    for t = 1:nT
        dm = pickMain(eM, sM, kM(t), P); dw = pickWav(eW, sW, kW(t), P);
        fprintf('   %5.2f/%-5.2f  ', sum(~bt(dm))/Tt, sum(~bt(dw))/Tt);
    end
    fprintf('   (main/db4, measured on unseen 45 s)\n');

    % ---- injection of transients of each amplitude, meter noise only ----
    [sh, pk] = shape(fs, tau1, tau2);
    rng(300 + f);
    fprintf('Se main/db4 at equal false-alarm budget (meter noise only, nothing added):\n');
    fprintf('  A p.u.  SNRi dB |'); fprintf('  FA<=%-5.2f    ', FAt); fprintf('\n');
    for a = 1:nA
        fprintf('  %5.2f  %6.1f  |', Amps(a), SNRtab(f,a));
        cnt = zeros(nT, 4);
        for rep = 1:nRep
            [xr, ev] = injectOnly(xt, bt, fs, Amps(a), sh, pk);
            [eM, sM] = scoreMain(xr, P);  [eW, sW] = scoreWav(xr);
            blk = floor((ev-1)/(Tblk*fs)) + 1;  nb = ceil(numel(xt)/(Tblk*fs));
            for t = 1:nT
                dm = pickMain(eM, sM, kM(t), P); dm = dm(~bt(dm));
                dw = pickWav(eW, sW, kW(t), P);  dw = dw(~bt(dw));
                hm = hitMask(dm, ev, P.tol);  hw = hitMask(dw, ev, P.tol);
                B = zeros(nb, 4);
                for j = 1:numel(ev)
                    B(blk(j), :) = B(blk(j), :) + [hm(j) ~hm(j) hw(j) ~hw(j)];
                end
                POOL{a,t} = [POOL{a,t}; B];
                cnt(t,:) = cnt(t,:) + sum(B, 1);
            end
        end
        for t = 1:nT
            fprintf('  %4.2f/%-4.2f    ', cnt(t,1)/sum(cnt(t,1:2)), cnt(t,3)/sum(cnt(t,3:4)));
        end
        fprintf('\n');
    end
end

fprintf('\n#####################################################################\n');
fprintf('POOLED over all files: Se(main) - Se(db4)  [95%% block-bootstrap CI]\n');
fprintf('  A p.u. |'); fprintf('  FA<=%-5.2f                 ', FAt); fprintf('\n');
for a = 1:nA
    fprintf('  %5.2f  |', Amps(a));
    for t = 1:nT
        B = POOL{a,t}; B = B(sum(B,2) > 0, :);
        if isempty(B), fprintf('  %-26s', 'n/a'); continue; end
        se = @(B) [sum(B(:,1))/sum(B(:,1)+B(:,2)), sum(B(:,3))/sum(B(:,3)+B(:,4))];
        s0 = se(B); d = zeros(2000,1); n = size(B,1);
        for b = 1:2000, s = se(B(randi(n,n,1),:)); d(b) = s(1)-s(2); end
        lo = myq(d, 0.025); hi = myq(d, 0.975);
        mk = '='; if lo > 0, mk = '+'; elseif hi < 0, mk = '-'; end
        fprintf('  %4.2f-%4.2f=%+.2f[%+.2f,%+.2f]%s', s0(1), s0(2), s0(1)-s0(2), lo, hi, mk);
    end
    fprintf('\n');
end
fprintf('(+ = D-TMA significantly more sensitive at the same false-alarm rate)\n');
fprintf('\nintrinsic SNR (dB) of each amplitude per file:\n');
for f = 1:nF, fprintf('  %-17s', files{f,2}); fprintf(' %6.1f', SNRtab(f,:)); fprintf('\n'); end
fprintf('\nTotal run time: %.1f s\n==== END ====\n', toc);

% =========================================================================
%                              FUNCTIONS
% =========================================================================
function [e, s] = scoreMain(x, P)
xp = dtma(x, P); mu = movmedian(xp, P.L); e = abs(xp - mu);
s  = max(1.4826*movmedian(e, P.L), eps);
end

function idx = pickMain(e, s, k, P)
score = e./s; idx = lmax(e, k*s);
idx = mergeg(idx, score, P.w1); idx = mergeg(idx, score, P.Tmin);
end

function [e, sg] = scoreWav(x)
g = [-0.0105974018 -0.0328830117 0.0308413818 0.1870348117 ...
     -0.0279837694 -0.6308807679 0.7148465706 -0.2303778133];
e = abs(conv(x(:), g(:), 'same')); sg = max(median(e)/0.6745, eps);
end

function idx = pickWav(e, sg, k, P)
idx = lmax(e, k*sg); idx = mergeg(idx, e, P.Tmin);
end

function h = hitMask(det, ev, tol)
det = det(:); used = false(size(det)); h = false(numel(ev), 1);
for j = 1:numel(ev)
    if isempty(det), break; end
    d = abs(det - ev(j)); d(used) = inf; [dm, k] = min(d);
    if dm <= tol, used(k) = true; h(j) = true; end
end
end

function [x, ev] = injectOnly(bg, bseg, fs, A0, sh, pk)
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
x = bg(:) + s;
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


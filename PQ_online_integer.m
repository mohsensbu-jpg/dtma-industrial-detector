% =========================================================================
% PQ_online_integer.m
% ONLINE (causal, sample-by-sample, integer) version of Mohsen Karimi's
% D-TMA + local rule, tested on the LOW-COST meter data (BidyutPi, 10 kHz).
%   filter : D-TMA-I, causal running sums only, integer input (ADC-like)
%            y = S_K2( K1^2 * q[n-(K1-1)] - S_K1(S_K1 q) )
%   noise  : local scale = recursive mean |y| (EMA, time constant ~10
%            cycles = 2^m samples, update = shift), winsorised at 4x
%            so events do not inflate it  (replaces the moving median,
%            which is too heavy for a microcontroller)
%   decide : local max of |y| above k*s, merge within G = max(w1,Tmin);
%            a detection is emitted G samples after its peak (online)
%   start  : 0.2 s warm-up, then 15 s calibration (nothing injected)
% Compared with the OFFLINE version (PQ_lowcost_equal_FA.m) at the same false-alarm budget.
% Reports: Se, false alarms/s, decision latency, ops/sample, memory,
% largest internal integer (overflow check).
% Text output only, no toolbox. Uses the 4 CSV files already in MATLAB Drive.
% =========================================================================
clear; clc; tic;
fprintf('==== PQ_online_integer  (%s) ====\n', datestr(now));
files = {'V_I_Waveforms_330W_Resistive_Load.csv',          'resistive 330 W';
         'V_I_Waveforms_0.5HP_Induction_Motor_No_Load.csv', 'induction motor';
         'V_I_Waveforms_50W_Ceiling_Fan_Aged.csv',          'ceiling fan 50 W';
         'V_I_Waveforms_SMPS_Mobile_Charger.csv',           'SMPS charger'};
Amps = [0.3 0.2 0.15];  FAt = [0.3 0.1];  Tcal = 15;  nRep = 3;
QS = 2048;  Gc = 4;                       % input scale (1 p.u. = 2048 counts), winsor factor
Wsec = 0.5e-3;  tau1 = 100e-6;  tau2 = 400e-6;
nA = numel(Amps); nT = numel(FAt); nF = size(files,1);
SUM = zeros(nA, nT, 4);  LAT = [];  MXall = 0;

for f = 1:nF
    fn = files{f,1};
    fprintf('\n#####################################################################\n');
    fprintf('FILE %d: %s\n', f, files{f,2});
    if ~isfile(fn), fprintf('file missing - run PQ_four_loads.m first\n'); continue; end
    [x, fs, f0, bad] = loadUniform(fn);
    N = numel(x);  P = params(fs, Wsec, f0);
    K1 = ceil(P.w1/2);  K2 = ceil(P.w2/2);  m = round(log2(P.L));
    D = (K1-1) + (K2-1);  G = max(P.w1, P.Tmin);  warm = round(0.2*fs);
    if f == 1
        fprintf('K1=%d K2=%d  filter delay D=%d smp (%.2f ms)  merge/emit G=%d smp (%.2f ms)  EMA 2^-%d (%.0f ms)\n', ...
            K1, K2, D, D/fs*1e3, G, G/fs*1e3, m, 2^m/fs*1e3);
    end
    q  = round(x*QS);
    ic = 1:round(Tcal*fs);  it = (round(Tcal*fs)+1):N;
    Tc = sum(~bad(ic))/fs;  ks = 1:0.05:30;

    % ---- calibration: online and offline, nothing injected ----
    [y, s, mx] = streamFilter(q(ic), K1, K2, m, Gc, warm);  MXall = max(MXall, mx);
    [eM, sM] = scoreMain(x(ic), P);
    kO = nan(1,nT); kF = nan(1,nT);
    for kk = ks
        if any(isnan(kO))
            d = decide(y, s, kk, G, warm) - D; d = d(d >= 1);
            fa = sum(~bad(ic(d)))/Tc;  kO(isnan(kO) & fa <= FAt) = kk;
        end
        if any(isnan(kF)) && kk >= 3
            fa = sum(~bad(ic(pickMain(eM, sM, kk, P))))/Tc;  kF(isnan(kF) & fa <= FAt) = kk;
        end
        if ~any(isnan([kO kF])), break; end
    end
    kO(isnan(kO)) = ks(end); kF(isnan(kF)) = ks(end);

    % ---- false alarms on unseen test part, nothing injected ----
    xt = x(it); qt = q(it); bt = bad(it); Tt = sum(~bt)/fs;
    [y, s, mx] = streamFilter(qt, K1, K2, m, Gc, warm);  MXall = max(MXall, mx);
    [eM, sM] = scoreMain(xt, P);
    fprintf('budget FA/s     :'); fprintf('   %5.2f          ', FAt); fprintf('\n');
    fprintf('k online/offline:'); fprintf('   %5.2f/%-5.2f    ', [kO; kF]); fprintf('\n');
    fprintf('test FA/s on/off:');
    for t = 1:nT
        d = decide(y, s, kO(t), G, warm) - D; d = d(d >= 1);
        dF = pickMain(eM, sM, kF(t), P);
        fprintf('   %5.2f/%-5.2f    ', sum(~bt(d))/Tt, sum(~bt(dF))/Tt);
    end
    fprintf('  (unseen 45 s)\n');

    % ---- injected transients, meter noise only ----
    [sh, pk] = shape(fs, tau1, tau2);  rng(300 + f);
    fprintf('Se online/offline at equal budget, and online decision latency:\n');
    for a = 1:nA
        cnt = zeros(nT, 4);  lat = [];
        for rep = 1:nRep
            [xr, ev] = injectOnly(xt, bt, fs, Amps(a), sh, pk);
            [y, s, mx] = streamFilter(round(xr*QS), K1, K2, m, Gc, warm);  MXall = max(MXall, mx);
            [eM, sM] = scoreMain(xr, P);
            for t = 1:nT
                [dO, tO] = decide(y, s, kO(t), G, warm);
                pos = dO - D;  keep = pos >= 1;  pos = pos(keep);  tO = tO(keep);
                keep = ~bt(pos);  pos = pos(keep);  tO = tO(keep);
                [hO, jO] = hitIdx(pos, ev, P.tol);
                dF = pickMain(eM, sM, kF(t), P);  dF = dF(~bt(dF));
                hF = hitIdx(dF, ev, P.tol);
                cnt(t,:) = cnt(t,:) + [sum(hO) sum(~hO) sum(hF) sum(~hF)];
                if t == nT, lat = [lat; tO(jO(hO)) - ev(hO)]; end %#ok<AGROW>
            end
        end
        SUM(a,:,:) = SUM(a,:,:) + reshape(cnt, [1 nT 4]);
        LAT = [LAT; lat]; %#ok<AGROW>
        fprintf('  A=%.2f p.u. |', Amps(a));
        for t = 1:nT
            fprintf('  FA<=%.2f: %.2f/%.2f', FAt(t), cnt(t,1)/sum(cnt(t,1:2)), cnt(t,3)/sum(cnt(t,3:4)));
        end
        if ~isempty(lat)
            fprintf('  | latency med/q95 = %.1f/%.1f ms', myq(lat,0.5)/fs*1e3, myq(lat,0.95)/fs*1e3);
        end
        fprintf('\n');
    end
end

fprintf('\n#####################################################################\nSUMMARY (all files pooled)\n');
for a = 1:nA
    fprintf('  A=%.2f p.u. |', Amps(a));
    for t = 1:nT
        c = squeeze(SUM(a,t,:));
        fprintf('  FA<=%.2f: Se online %.3f  offline %.3f |', FAt(t), c(1)/(c(1)+c(2)), c(3)/(c(3)+c(4)));
    end
    fprintf('\n');
end
if ~isempty(LAT)
    fprintf('decision latency (event -> alarm), online: median %.1f ms, q95 %.1f ms, max %.1f ms\n', ...
        myq(LAT,0.5)/fs*1e3, myq(LAT,0.95)/fs*1e3, max(LAT)/fs*1e3);
end
fprintf('largest |internal integer| seen: %d  (int32 limit 2147483647 -> %s)\n', ...
    MXall, ternary(MXall < 2^31-1, 'OK', 'OVERFLOW'));
fprintf('cost per sample (K1=%d, K2=%d): 13 add/sub, 2 multiplies (K1^2*q, k*s), 2 shifts, 5 compares\n', K1, K2);
fprintf('memory: %d integers of buffers (2*K1 + 2*K2 + K1-1) + 10 state variables\n', 2*K1 + 2*K2 + K1 - 1);
fprintf('\nTotal run time: %.1f s\n==== END ====\n', toc);

% =========================================================================
%                              FUNCTIONS
% =========================================================================
function [y, s, mx] = streamFilter(q, K1, K2, m, Gc, warm)
% genuine sample-by-sample causal loop, integer values only
N = numel(q); y = zeros(N,1); s = ones(N,1);
b1 = zeros(K1,1); c1 = zeros(K1,1); b2 = zeros(K2,1); c2 = zeros(K2,1);
Dl = K1 - 1; dl = zeros(max(Dl,1),1);
S1 = 0; S2 = 0; T1 = 0; T2 = 0; i1 = 1; i2 = 1; id = 1;
acc = 0; sumA = 0; mx = 0; M2 = 2^m;
for n = 1:N
    xn = q(n);
    S1 = S1 + xn - b1(i1);  b1(i1) = xn;          % running sum 1
    S2 = S2 + S1 - c1(i1);  c1(i1) = S1;          % running sum 2 (TMA*K1^2)
    i1 = i1 + 1; if i1 > K1, i1 = 1; end
    if Dl > 0, xd = dl(id); dl(id) = xn; id = id + 1; if id > Dl, id = 1; end
    else, xd = xn; end
    h  = K1*K1*xd - S2;                           % exact DC rejection
    T1 = T1 + h - b2(i2);  b2(i2) = h;
    T2 = T2 + T1 - c2(i2); c2(i2) = T1;
    i2 = i2 + 1; if i2 > K2, i2 = 1; end
    e = abs(T2);
    if abs(S2) > mx, mx = abs(S2); end
    if abs(h)  > mx, mx = abs(h);  end
    if e       > mx, mx = e;       end
    if n <= warm                                  % warm-up: plain mean |y|
        sumA = sumA + e;  acc = floor(sumA*M2/n);
    else
        sc  = floor(acc/M2);                      % = acc >> m
        acc = acc + min(e, Gc*sc) - sc;           % winsorised EMA
    end
    if acc > mx, mx = acc; end
    y(n) = T2;  s(n) = max(floor(acc/M2), 1);
end
end

function [det, tem] = decide(y, s, k, G, warm)
% same result as the per-sample loop: local max above k*s, chain merge
% within G samples (keep larger e/s), emit when G samples pass without a
% stronger peak. tem = emission sample (decision time).
e = abs(y); n = numel(e); i = (2:n-1)';
c = i(e(i) > e(i-1) & e(i) >= e(i+1) & e(i) > k*s(i) & i > warm);
det = zeros(numel(c),1); tem = det; nd = 0; pend = 0; ps = 0;
for j = 1:numel(c)
    i = c(j); sc = e(i)/s(i);
    if pend > 0 && i - pend < G
        if sc > ps, pend = i; ps = sc; end
    else
        if pend > 0, nd = nd + 1; det(nd) = pend; tem(nd) = pend + G; end
        pend = i; ps = sc;
    end
end
if pend > 0, nd = nd + 1; det(nd) = pend; tem(nd) = pend + G; end
det = det(1:nd); tem = tem(1:nd);
end

function [h, jdx] = hitIdx(det, ev, tol)
det = det(:); used = false(size(det)); h = false(numel(ev),1); jdx = ones(numel(ev),1);
for j = 1:numel(ev)
    if isempty(det), break; end
    d = abs(det - ev(j)); d(used) = inf; [dm, k] = min(d);
    if dm <= tol, used(k) = true; h(j) = true; jdx(j) = k; end
end
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end

function [e, s] = scoreMain(x, P)
xp = dtma(x, P); mu = movmedian(xp, P.L); e = abs(xp - mu);
s  = max(1.4826*movmedian(e, P.L), eps);
end

function idx = pickMain(e, s, k, P)
score = e./s; idx = lmax(e, k*s);
idx = mergeg(idx, score, P.w1); idx = mergeg(idx, score, P.Tmin);
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

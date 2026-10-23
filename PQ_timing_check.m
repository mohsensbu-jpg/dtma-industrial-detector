% =========================================================================
% PQ_timing_check.m   (run PQ_four_loads.m first (it downloads the data))
% Question: is the high "real noise" (intrinsic SNR 7 dB) and the
% 16 detections/s on a clean resistive load caused by NON-UNIFORM
% sampling of the Raspberry-Pi meter, not by the grid?
% 1) statistics of the sampling interval dt
% 2) are the strongest detections at sampling gaps?
% 3) resample onto a uniform grid using the real timestamps and repeat:
%    residual noise and detections per second (before / after)
% Text output only, no toolbox.
% =========================================================================
clear; clc;
fn = 'V_I_Waveforms_330W_Resistive_Load.csv';
T  = readtable(fn, 'VariableNamingRule', 'preserve');
t  = double(T.time_s(:));  vV = double(T.voltage_V(:));  vA = double(T.adc_voltage(:));
N  = numel(t);
fprintf('==== PQ_timing_check ====\nN = %d   t(1) = %.6f   t(end) = %.6f s\n', N, t(1), t(end));

%% 1) sampling interval
dt = diff(t);
fprintf('\n-- 1) sampling interval dt (us) --\n');
fprintf('non-increasing steps (dt<=0): %d\n', sum(dt <= 0));
p = [0 0.001 0.01 0.05 0.25 0.5 0.75 0.95 0.99 0.999 1];
for k = 1:numel(p), fprintf('  q%-6.3f = %9.1f us\n', p(k), myq(dt, p(k))*1e6); end
fprintf('mean dt = %.2f us  -> mean fs = %.1f Hz\n', mean(dt)*1e6, 1/mean(dt));
fprintf('dt > 2x median: %d  (%.3f %%)   dt > 5x median: %d   dt < 0.5x median: %d\n', ...
    sum(dt > 2*median(dt)), 100*mean(dt > 2*median(dt)), sum(dt > 5*median(dt)), sum(dt < 0.5*median(dt)));
[~, o] = sort(dt, 'descend');
fprintf('largest gaps: '); fprintf('%.4f s (%.0f us)  ', [t(o(1:5))'; dt(o(1:5))'*1e6]); fprintf('\n');
fprintf('voltage_V vs adc_voltage linear? corr = %.6f\n', corrOwn(vV, vA));

%% 2) analysis on the sample index (no timing correction) and on a uniform grid
fsU = round(1/median(dt));             % nominal rate
tu  = (t(1) : 1/fsU : t(end))';
vu  = interp1(t, vV, tu, 'linear');
cases = {'sample index (no timing fix), voltage_V', vV, fsU;
         'UNIFORM grid from timestamps',        vu, fsU};
for c = 1:2
    x = cases{c,2}; x = x - mean(x); fs = cases{c,3};
    [tc, nUp] = upcross(x); f0 = (nUp-1)/(tc(end)-tc(1))*fs;
    x = x / (sqrt(2)*sqrt(mean(x.^2)));
    nz = harmonicResidual(x, fs, tc);
    P  = params(fs, 0.5e-3, f0);
    [d, sc] = detMain(x, P);
    fprintf('\n-- 2.%d) %s --\n', c, cases{c,1});
    fprintf('N = %d  f0 = %.3f Hz  residual noise RMS = %.5f p.u. (intrinsic SNR for 0.2 p.u. = %.1f dB)\n', ...
        numel(x), f0, sqrt(mean(nz.^2)), 20*log10(0.2/sqrt(mean(nz.^2))));
    fprintf('main detector: %d detections (%.2f /s)\n', numel(d), numel(d)/(numel(x)/fs));
    if c == 1
        % distance of each detection to the nearest large sampling gap
        g = find(dt > 2*median(dt));           % gap after sample g
        if ~isempty(g) && ~isempty(d)
            dist = zeros(numel(d),1);
            for j = 1:numel(d), dist(j) = min(abs(g - d(j))); end
            fprintf('detections within 3 samples of a gap (dt>2x median): %d of %d (%.1f %%)\n', ...
                sum(dist <= 3), numel(d), 100*mean(dist <= 3));
            fprintf('fraction of all samples within 3 samples of a gap: %.1f %% (chance level)\n', ...
                100*mean(minDist((1:N)', g) <= 3));
        end
        [~, o] = sort(sc(d), 'descend');
        fprintf('strongest 5: ');
        for j = o(1:min(5,end))', fprintf('[t=%.4f e/s=%.1f dt_near=%.0fus] ', t(d(j)), sc(d(j)), ...
                1e6*max(dt(max(1,d(j)-3):min(N-1,d(j)+3)))); end
        fprintf('\n');
    end
end
fprintf('\n==== END ====\n');

% ---------------- functions (shared) -----------------------------
function md = minDist(a, g)
g = sort(g(:)); md = inf(size(a));
k = discretize(a, [-inf; g; inf]);       % interval index
for i = 1:numel(a)
    if k(i) > 1, md(i) = min(md(i), a(i) - g(k(i)-1)); end
    if k(i) <= numel(g), md(i) = min(md(i), g(k(i)) - a(i)); end
end
end
function c = corrOwn(a, b)
a = a - mean(a); b = b - mean(b); c = sum(a.*b)/sqrt(sum(a.^2)*sum(b.^2));
end
function P = params(fs, Wsec, f0)
P.fs = fs; P.W = max(3, round(Wsec*fs)); P.w1 = max(3, round(3.3*P.W));
P.w2 = max(3, round(0.06*P.W)); P.L = round(10/f0*fs); P.Tmin = round(2e-3*fs);
end
function [idx, score] = detMain(x, P)
xp = tma(x(:) - tma(x(:), P.w1), P.w2);
mu = movmedian(xp, P.L); e = abs(xp - mu);
s  = max(1.4826*movmedian(e, P.L), eps); score = e./s;
idx = lmax(e, sqrt(2*log(P.L))*s);
idx = mergeg(idx, score, P.w1); idx = mergeg(idx, score, P.Tmin);
end
function y = tma(x, w)
K = ceil(w/2); y = movmean(x, K);
if mod(K,2) == 0, y = movmean(y, [K/2-1 K/2]); else, y = movmean(y, K); end
end
function idx = lmax(e, thr)
n = numel(e); i = (2:n-1)';
idx = i(e(i) > e(i-1) & e(i) >= e(i+1) & e(i) > thr(i));
end
function out = mergeg(idx, score, gap)
out = zeros(numel(idx),1); c = 0;
for j = 1:numel(idx)
    i = idx(j);
    if c > 0 && i - out(c) < gap
        if score(i) > score(out(c)), out(c) = i; end
    else, c = c + 1; out(c) = i; end
end
out = out(1:c);
end
function q = myq(x, p)
x = sort(x(:)); n = numel(x); pos = 1 + (n-1)*p; lo = floor(pos); hi = ceil(pos);
q = x(lo) + (pos-lo)*(x(hi) - x(lo));
end
function [tc, nUp] = upcross(x)
a = 0.1*max(abs(x)); st = 0; tc = zeros(ceil(numel(x)/2),1); nUp = 0; last = 2;
for i = 2:numel(x)
    if x(i-1) <= 0 && x(i) > 0, last = i; end
    if st <= 0 && x(i) > a
        if st < 0, nUp = nUp + 1; k = last; tc(nUp) = (k-1) + (-x(k-1))/(x(k)-x(k-1)); end
        st = 1;
    elseif st >= 0 && x(i) < -a, st = -1; end
end
tc = tc(1:nUp);
end
function nz = harmonicResidual(v, fs, tc)
N = numel(v); B = round(0.2*fs); nz = zeros(N,1);
f0g = (numel(tc)-1)/(tc(end)-tc(1))*fs;
for b = 1:floor(N/B)
    ii = ((b-1)*B+1 : b*B)'; c = tc(tc >= ii(1)-1 & tc <= ii(end)-1);
    if numel(c) >= 3, f0 = (numel(c)-1)/(c(end)-c(1))*fs; else, f0 = f0g; end
    tt = (0:B-1)'/fs; nh = min(floor(0.45*fs/f0), 100); H = ones(B, 1+2*nh);
    for h = 1:nh, H(:,2*h) = cos(2*pi*h*f0*tt); H(:,2*h+1) = sin(2*pi*h*f0*tt); end
    nz(ii) = v(ii) - H*(H\v(ii));
end
nz = nz(1:floor(N/B)*B);
end

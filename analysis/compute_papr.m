function papr_dB = compute_papr(x)
% 根据课上讲的PAPR的定义计算信号的PAPR（dB）
    peak_power = max(abs(x).^2);
    avg_power  = mean(abs(x).^2);
    papr_dB    = 10 * log10(peak_power / avg_power);
end

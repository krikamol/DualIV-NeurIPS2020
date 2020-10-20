function y = storeg(x_latent, price, psd, pmu, ymu, ysd)
% emocell = num2cell(emotion_feature, 2);
% emoc = cellfun(@emocoef, emocell);
time = x_latent(:, 1);
emoc = x_latent(:, 2);
g = sensf(time) .* emoc * 10. + (emoc .* sensf(time) - 2.0) .* (psd * price + pmu);
y = (g - ymu) / ysd;
end
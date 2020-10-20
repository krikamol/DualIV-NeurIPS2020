function emoc = emocoef(emo)
% emo: the one-hot representation of the emotion index
% emoc: the numerical representation of the emotion
emoc = dot(emo, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]);
end
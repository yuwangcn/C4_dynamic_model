function LeafMs = C4leafMetaIni()

leafInis = zeros(1,7);
leafInis= LeafIni;
MetaInis = zeros(1,87);
MetaInis= C4Ini;

for m = 1:7
    LeafMs(m) = leafInis(m);
end

for m=1:87
    LeafMs(7+m)= MetaInis(m);
end
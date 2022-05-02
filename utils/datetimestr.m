function s = datetimestr()
%DATETIMESTR returns the date and time in string format.

dati = datetime;
s = sprintf("%d/%02d/%02d-%02d:%02d:%02d", dati.Year, dati.Month, dati.Day, ...
                                           dati.Hour, dati.Minute, round(dati.Second));
end
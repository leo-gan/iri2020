import iri2020
time = "2015-12-13T10"
altkmrange = (100, 1000, 10.0)
glat = 65.1
glon = -147.5
iri = iri2020.IRI(time, altkmrange, glat, glon)
print("altkm length:", len(iri["alt_km"]))
print("ne length:", len(iri["ne"]))
print("ne values (first 20):", iri["ne"][:20].values)
print("NmF2:", iri.NmF2.item())
print("hmF2:", iri.hmF2.item())

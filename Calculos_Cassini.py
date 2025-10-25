import numpy as np


# Constantes
#* AU^3/dia^2
μ_Sol = 0.0002959122

#* Constante μ [AU^3/dia^2]
μ_planetas = {
    "Terra": 3.986e-14,
    "Venus": 3.24859e-14,
    "Jupiter": 1.26686534e-13
}

#* AU
R_Orb_planetas = {
    "Venus": 0.723,
    "Terra": 1,
    "Jupiter": 5.204
}


Vetores_Cassini = {}

with open("Cassini_Vetor.csv", "r") as f:
    for line in f:
        if line.startswith("Vetor_"):
            split_line = line.strip().split(";")
            secao = split_line[0].replace("Vetor_", "").strip()
            encontroucabecalho = False

        elif line.startswith("X"):
            cabecalho = split_line[0].replace("X", "").strip()
            encontroucabecalho = True

        elif encontroucabecalho:
            valores = float(split_line[1].replace(",",".").strip())
            vetores = []
            for i in range (6):
                vetor = float(line.split(";")[i].strip().replace(",","."))
                vetores.append(vetor)
            rvec = np.array(vetores[0:3])
            vvec = np.array(vetores[3:6])

# Objetos Auxiliares
class Orbita:
    def vet_estado(self, μ_Sol, rvec, vvec):
        # Magnitude Posição e Velocidade
        r = np.linalg.norm(rvec)
        v = np.linalg.norm(vvec)

        # Semieixo maior
        self.a = 1 / ((2 / r) - (v**2 / μ_Sol))

        # Vetor momento angular
        hvec = np.cross(rvec, vvec)

        # Magnitude momento angular
        h = np.linalg.norm(hvec)

        # Semilado reto e Excentricidade
        self.p = h**2 / μ_Sol
        self.e = np.sqrt(1 - (self.p / self.a))


    def raio_orbita(self, f):
        return (self.p) / (1 + self.e * np.cos(f))

    def pos_orb(self, ω, f):
        xpq = self.raio_orbita(f) * np.cos(ω + f)
        ypq = self.raio_orbita(f) * np.sin(ω + f)
        return xpq, ypq

    def v_orb(self, μ, f):
        vxpq = -np.sqrt((μ / self.p)) * np.sin(f)
        vypq = np.sqrt((μ / self.p)) * (self.e + np.cos(f))
        return vxpq, vypq





import numpy as np
import matplotlib.pyplot as plt

μ_Sol = 0.0002959122

μ_planetas = {
    "Terra": 3.986e-14,
    "Venus": 3.24859e-14,
    "Jupiter": 1.26686534e-13
}

class Orbita:
    def __init__(self):
        self.secao_nome = ""
        self.a = 0
        self.e = 0
        self.p = 0
        self.omega = 0
        self.rvec = np.zeros(3)
        self.vvec = np.zeros(3)
        self.e_vec = np.zeros(3)

    def vet_estado(self, μ_Sol, rvec, vvec):
        self.rvec = rvec
        self.vvec = vvec

        r = np.linalg.norm(rvec)
        v = np.linalg.norm(vvec)

        self.a = 1 / ((2 / r) - (v**2 / μ_Sol))

        hvec = np.cross(rvec, vvec)

        termo1 = (v**2 - (μ_Sol / r)) * rvec
        termo2 = np.dot(rvec, vvec) * vvec
        self.e_vec = (1 / μ_Sol) * (termo1 - termo2)

        self.e = np.linalg.norm(self.e_vec)
        self.p = np.linalg.norm(hvec)**2 / μ_Sol

        self.omega = np.arctan2(self.e_vec[1], self.e_vec[0])

    def raio_orbita(self, f):
        return (self.p) / (1 + self.e * np.cos(f))

    def pos_orb(self, ω, f):
        xpq = self.raio_orbita(f) * np.cos(ω + f)
        ypq = self.raio_orbita(f) * np.sin(ω + f)
        return xpq, ypq

    def get_true_anomaly(self, r_target):
        e_unit = self.e_vec / self.e
        h_vec = np.cross(self.rvec, self.vvec)
        h_unit = h_vec / np.linalg.norm(h_vec)
        q_unit = np.cross(h_unit, e_unit)

        x_proj = np.dot(r_target, e_unit)
        y_proj = np.dot(r_target, q_unit)

        return np.arctan2(y_proj, x_proj)

Vetores_Cassini = {}
Vetores_Planetas = {}

encontroucabecalho = False
secao = None

try:
    with open("Cassini_Vetor.csv", "r") as f:
        for line in f:
            line_branco = line.strip()
            if not line_branco or line_branco.startswith(";") or line_branco.startswith(" "):
                continue
            if line.startswith("Vetor_"):
                split_line = line.strip().split(";")
                secao = split_line[0].strip()
                encontroucabecalho = False
            elif line.startswith("X"):
                encontroucabecalho = True
            elif encontroucabecalho and secao:
                vetores = []
                try:
                    for i in range (6):
                        vetor = float(line.split(";")[i].strip().replace(",","."))
                        vetores.append(vetor)
                    rvec = np.array(vetores[0:3])
                    vvec = np.array(vetores[3:6])

                    nova_orbita = Orbita()
                    nova_orbita.secao_nome = secao
                    nova_orbita.vet_estado(μ_Sol, rvec, vvec)
                    Vetores_Cassini[secao] = nova_orbita
                    encontroucabecalho = False
                except Exception:
                    encontroucabecalho = False
except FileNotFoundError:
    print("O arquivo Cassini_Vetor.csv não foi encontrado.")

try:
    with open("Vetores_Planetas.csv", "r") as f:
        for line in f:
            line_branco = line.strip()
            if not line_branco or line_branco.startswith(";") or line_branco.startswith(" "):
                continue
            if line.startswith("Vetores_"):
                split_line = line.strip().split(";")
                secao = split_line[0].strip()
                encontroucabecalho = False
            elif line.startswith("X"):
                encontroucabecalho = True
            elif encontroucabecalho and secao:
                vetores = []
                try:
                    for i in range (6):
                        vetor = float(line.split(";")[i].strip().replace(",","."))
                        vetores.append(vetor)
                    rvec = np.array(vetores[0:3])
                    vvec = np.array(vetores[3:6])

                    nova_orbita = Orbita()
                    nova_orbita.secao_nome = secao
                    nova_orbita.vet_estado(μ_Sol, rvec, vvec)
                    Vetores_Planetas[secao] = nova_orbita
                    encontroucabecalho = False
                except Exception:
                    encontroucabecalho = False
except FileNotFoundError:
    print("O arquivo Vetores_Planetas.csv não foi encontrado.")

plt.figure(figsize=(12, 10))
plt.title("Trajetória da Missão Cassini: Segmentos Sólidos e Órbitas Tracejadas")
plt.xlabel("x (AU)")
plt.ylabel("y (AU)")
plt.grid(True, linestyle='--', alpha=0.3)

plt.plot(0, 0, 'yo', markersize=12, label='Sol')

f_full = np.linspace(0, 2 * np.pi, 400)
for nome, orbita in Vetores_Planetas.items():
    x_vals, y_vals = orbita.pos_orb(orbita.omega, f_full)
    plt.plot(x_vals, y_vals, '--', color='gray', alpha=0.4, linewidth=0.8)
    plt.plot(orbita.rvec[0], orbita.rvec[1], 'o', markersize=6, label=nome)

cores = ['blue', 'green', 'red', 'orange', 'purple', 'cyan', 'magenta']
chaves_cassini = list(Vetores_Cassini.keys())

for i, nome in enumerate(chaves_cassini):
    orbita_atual = Vetores_Cassini[nome]
    cor_atual = cores[i % len(cores)]

    x_full, y_full = orbita_atual.pos_orb(orbita_atual.omega, f_full)
    plt.plot(x_full, y_full, ':', color=cor_atual, alpha=0.6, linewidth=1)

    f_start = orbita_atual.get_true_anomaly(orbita_atual.rvec)

    if i < len(chaves_cassini) - 1:
        orbita_proxima = Vetores_Cassini[chaves_cassini[i+1]]
        f_end = orbita_atual.get_true_anomaly(orbita_proxima.rvec)

        if f_end < f_start:
            f_end += 2 * np.pi
    else:
        f_end = f_start + np.pi

    f_segmento = np.linspace(f_start, f_end, 200)
    x_seg, y_seg = orbita_atual.pos_orb(orbita_atual.omega, f_segmento)

    plt.plot(x_seg, y_seg, '-', color=cor_atual, linewidth=2.5, label=nome)

    plt.plot(x_seg[0], y_seg[0], 'o', color=cor_atual, markersize=6, markeredgecolor='black')
    plt.plot(x_seg[-1], y_seg[-1], 's', color=cor_atual, markersize=6, markeredgecolor='black')

plt.axis('equal')
plt.legend(loc='upper right', fontsize='small', framealpha=0.9)
plt.tight_layout()
plt.show()

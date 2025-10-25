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

# Objetos Auxiliares
class Orbita:
    def __init__(self):
        self.secao_nome = ""
        self.a = 0
        self.e = 0
        self.p = 0
        self.rvec = np.zeros(3)
        self.vvec = np.zeros(3)

    def vet_estado(self, μ_Sol, rvec, vvec):
        self.rvec = rvec
        self.vvec = vvec

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

# Vetores da Cassini e Planetas
Vetores_Cassini = {}
Vetores_Planetas = {}

# Condições iniciais para leitura
encontroucabecalho = False
secao = None


try:
    with open("Cassini_Vetor.csv", "r") as f:

        for line in f:
            line_branco = line.strip()

            if not line_branco or line_branco.startswith(";") or line_branco.startswith(" "):
                continue
            # Separa os dados
            if line.startswith("Vetor_"):
                split_line = line.strip().split(";")
                secao = split_line[0].strip()
                encontroucabecalho = False

            elif line.startswith("X"):
                encontroucabecalho = True

            # Organiza vetores posição e velocidade

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
                except Exception as e:
                    print(f"Erro ao processar a linha: {line}. Erro: {e}")
                    encontroucabecalho = False
except FileNotFoundError:
    print("O arquivo não foi encontrado.")


try:
    with open("Vetores_Planetas.csv", "r") as f:

        for line in f:
            line_branco = line.strip()

            if not line_branco or line_branco.startswith(";") or line_branco.startswith(" "):
                continue
            # Separa os dados
            if line.startswith("Vetores_"):
                split_line = line.strip().split(";")
                secao = split_line[0].strip()
                encontroucabecalho = False

            elif line.startswith("X"):
                encontroucabecalho = True

            # Organiza vetores posição e velocidade

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
                except Exception as e:
                    print(f"Erro ao processar a linha: {line}. Erro: {e}")
                    encontroucabecalho = False
except FileNotFoundError:
    print("O arquivo não foi encontrado.")

# Teste
try:
    print(f"  Semi-eixo maior (a): {Vetores_Cassini['Vetor_Terra_Venus'].a:.4f} AU")
except KeyError:
    print(" A chave não foi encontrada no dicionário.")
try:
    print(f"  Semi-eixo maior (a): {Vetores_Planetas['Vetores_Terra'].a:.4f} AU")
except KeyError:
    print(" A chave não foi encontrada no dicionário.")



with open("Dados_Orbita_Cassini.txt", "w") as f:
    for secao_nome, orbita_obj in Vetores_Cassini.items():
        f.write(f"Dados para {secao_nome}:\n")
        f.write(f"  Semi-eixo maior (a): {orbita_obj.a:.4f} AU\n")
        f.write(f"  Excentricidade (e): {orbita_obj.e:.4f}\n")
        f.write(f"  Semilado reto (p): {orbita_obj.p:.4f} AU\n")
        f.write(f"  Posicao (r): {orbita_obj.rvec} AU\n")
        f.write(f"  Velocidade (v): {orbita_obj.vvec} AU/dia\n")
        f.write("\n")

with open("Dados_Orbita_Planetas.txt", "w") as f:
    for secao_nome, orbita_obj in Vetores_Planetas.items():
        f.write(f"Dados para {secao_nome}:\n")
        f.write(f"  Semi-eixo maior (a): {orbita_obj.a:.4f} AU\n")
        f.write(f"  Excentricidade (e): {orbita_obj.e:.4f}\n")
        f.write(f"  Semilado reto (p): {orbita_obj.p:.4f} AU\n")
        f.write(f"  Posicao (r): {orbita_obj.rvec} AU\n")
        f.write(f"  Velocidade (v): {orbita_obj.vvec} AU/dia\n")
        f.write("\n")



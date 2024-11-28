#===============================================================================
# Packages
#===============================================================================

# Liste des packages nécessaires
packages <- c("rgl", "pracma", "randtoolbox")

# Installation des packages manquants
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Chargement des packages
lapply(packages, library, character.only = TRUE)

#===============================================================================
# Estimation de l'intégrale de exp(sin(x * y))
#===============================================================================

# Définition des limites des variables
a <- 0    # xmin
b <- 0.5  # xmax
c <- 2    # ymin
d <- 3    # ymax

# Nombre de simulations
n_simulations <- 1000  # Total de points à simuler
n_iterations <- 500      # Nombre d'itérations pour chaque méthode

# Fonction à intégrer
f_xy <- function(x, y) {
  exp(sin(x * y))  # Calcul de f(x, y)
}

# Création de la grille de points pour le graphique 3D
valeurs_x <- seq(a, b, length.out = 100)
valeurs_y <- seq(c, d, length.out = 100)

# Calcul des valeurs Z en utilisant outer pour appliquer f_xy sur chaque combinaison de x et y
valeurs_z <- outer(valeurs_x, valeurs_y, f_xy)

# Affichage du graphique 3D
persp3d(valeurs_x, valeurs_y, valeurs_z,
        col = "lightblue",                # Couleur de la surface
        xlab = "X",                       # Étiquette de l'axe X
        ylab = "Y",                       # Étiquette de l'axe Y
        zlab = "f(x, y)",                # Étiquette de l'axe Z
        phi = 30,                         # Angle de la vue en élévation
        theta = 30,                       # Angle de la vue en rotation
        lit = TRUE,                       # Éclairage de la surface
        axes = TRUE)                     # Afficher les axes

#===============================================================================
# Méthode 1 : Méthode de volume
#===============================================================================

# Définition des limites pour la hauteur (z)
z_min <- 0  # Valeur minimale pour z
z_max <- exp(sin(b * d))  # Valeur maximale pour z

# Tirage de points dans l'espace de volume
z_ech <- matrix(runif(n_simulations * n_iterations, z_min, z_max), ncol = n_iterations)

# Tirage de points pour x et y
x_ech <- matrix(runif(n_simulations * n_iterations, a, b), ncol = n_iterations)
y_ech <- matrix(runif(n_simulations * n_iterations, c, d), ncol = n_iterations)

# Calculer f(x, y) pour les échantillons tirés
z_surface <- f_xy(x_ech, y_ech)

# Vérifier combien de points sont sous la surface
points_en_dessous_surface <- z_ech <= z_surface

# Approximation de l'intégrale par la méthode du volume
facteur_volume <- (b - a) * (d - c) * z_max
Tn_evol <- cumsum(points_en_dessous_surface[, 1]) / (1:n_simulations) * facteur_volume

# Calcul de Tn pour toutes les simulations
calculer_estimation_volume <- function(points) {
  mean(points) * facteur_volume
}
Tn_vecteur <- apply(points_en_dessous_surface, MARGIN = 2, FUN = calculer_estimation_volume)

#===============================================================================
# Méthode 2 : Méthode de Monte Carlo basique
#===============================================================================

# Calcul de f(x, y) pour chaque échantillon tiré
valeurs_z_mc <- f_xy(x_ech, y_ech)

# Approximation de l'intégrale par Monte Carlo basique
aire_rectangle <- (b - a) * (d - c)
Sn_evol <- cumsum(valeurs_z_mc[, 1]) / (1:n_simulations) * aire_rectangle

# Calcul de Sn pour toutes les simulations
calculer_estimation_mc <- function(valeurs) {
  mean(valeurs) * aire_rectangle
}
Sn_vecteur <- apply(valeurs_z_mc, MARGIN = 2, FUN = calculer_estimation_mc)

#===============================================================================
# Méthode 3 : Méthode de variance réduite
#===============================================================================

# Définir plusieurs valeurs pour y : c (correspondant à y_min), la médiane, et d (correspondant à y_max) 
y_fixe <- c(c, (c + d) / 2, d)  # Création d'un vecteur y_fixe contenant c (min), la médiane, et d (max)

# Calculer les valeurs de f(x, y_fixe) pour chaque combinaison de x et y
f_x_yfixe <- sapply(y_fixe, function(y) sapply(valeurs_x, function(x) f_xy(x, y)))  # Création d'une matrice où chaque colonne correspond à f(x, y) pour une valeur de y fixée

# Tracer la courbe pour y = c (correspondant à y_min) avec une ligne noire pointillée
plot(valeurs_x, f_x_yfixe[, 1], type = "l", col = "black", lwd = 2, lty = 2,
     xlab = "X", ylab = "f(X, Y)", 
     main = "Graphique de f(x, y) pour y : c, médiane, d",  # Titre plus clair
     ylim = c(1, 3),
     cex.main = 0.8,   # Taille du titre réduite
     cex.lab = 0.8,    # Taille des étiquettes des axes réduite
     cex.axis = 0.8)  # Définir les limites de l'axe Y de 1 à 3

# Ajouter la courbe pour y = d (correspondant à y_max) avec une ligne noire pointillée
lines(valeurs_x, f_x_yfixe[, 3], col = "black", lwd = 2, lty = 2)

# Colorier la région entre c (y_min) et d (y_max) en bleu transparent
polygon(c(valeurs_x, rev(valeurs_x)), 
        c(f_x_yfixe[, 1], rev(f_x_yfixe[, 3])), 
        col = rgb(0.1, 0.1, 0.8, 0.2), border = NA)  # Remplissage en bleu transparent

# Ajouter la courbe pour la valeur médiane de y avec une ligne rouge
lines(valeurs_x, f_x_yfixe[, 2], col = "red", lwd = 2)

# Calculer la fonction p1 à tracer en bleu
p1_valeurs <- 8 / 7 * (1 + 3 * valeurs_x)  # Définir p1 en fonction de x, selon la formule donnée
lines(valeurs_x, p1_valeurs, col = "blue", lwd = 2)  # Tracer la ligne de p1 en bleu

# Ajouter une légende pour identifier les courbes
legend("bottomright", 
       legend = c(paste("y =", c),  # Indiquer la valeur de c (y_min)
                  paste("y =", (c + d) / 2),  # Indiquer la valeur médiane
                  paste("y =", d),  # Indiquer la valeur de d (y_max)
                  "p1 = 8/7(1 + 3x)"),  # Indiquer la formule de p1
       col = c("black", "red", "black", "blue"),  # Couleurs pour la légende
       lwd = 2, lty = c(2, 1, 2, 1))  # Styles de lignes pour la légende
grid()


# Définir p1, la densité de probabilité pour X
p1 <- function(x) 8 / 7 * (1 + 3 * x)

# Fonction de répartition de p1 (cumulée)
F1 <- function(t) {
  return(4 / 7 * (2 * t + 3 * t^2))  # Fonction F1(x) correspondant à la fonction de répartition de p1
}

# Inverse de la fonction de répartition F1 (permet de générer des échantillons de X selon p1)
F1_inverse <- function(z) {
  return((-1 + sqrt(1 + 21 / 4 * z)) / 3)  # Résolution analytique de l'inverse de F1
}

# Histogramme des valeurs générées par la fonction inverse de F1 avec une distribution uniforme
hist(F1_inverse(runif(10000)), 
     col = rgb(0.1, 0.5, 0.8, 0.7),  # Couleur des barres avec transparence
     border = "white",  # Couleur de la bordure des barres
     xlab = "Valeurs de X",  # Étiquette de l'axe X
     ylab = "Fréquence",  # Étiquette de l'axe Y
     main = paste("Histogramme des échantillons de X", "\nselon la densité de proba 1"),  # Titre avec saut de ligne
     xlim = c(0, 0.5),  # Limites de l'axe X (ajuster selon votre intervalle)
     ylim = c(0, 1400),  # Limites de l'axe Y (ajuster selon votre distribution)
     las = 1,  # Orientation des étiquettes de l'axe Y
     freq = TRUE,
     cex.main = 0.8,   # Taille du titre réduite
     cex.lab = 0.8,    # Taille des étiquettes des axes réduite
     cex.axis = 0.8)  # Fréquence ou densité (TRUE pour fréquence)


# Définir plusieurs valeurs pour x : a (correspondant à x_min), la médiane, et b (correspondant à x_max)
x_fixe <- c(a, (a + b) / 2, b)  # Création d'un vecteur x_fixe contenant a (min), la médiane, et b (max)

# Calculer les valeurs de f(x_fixe, y) pour chaque combinaison de y et x
f_xfixe_y <- sapply(x_fixe, function(x) sapply(valeurs_y, function(y) f_xy(x, y)))  # Création d'une matrice où chaque colonne correspond à f(x, y) pour une valeur de x fixée

# Tracer la courbe pour x = a (correspondant à x_min) avec une ligne noire pointillée
plot(valeurs_y, f_xfixe_y[, 1], type = "l", col = "black", lwd = 2, lty = 2,
     xlab = "Y", ylab = "f(X, Y)", 
     main = "Graphique de f(x, y) pour x : a, médiane, b",  # Titre plus clair
     ylim = c(1, 3),
     cex.main = 0.8,   # Taille du titre réduite
     cex.lab = 0.8,    # Taille des étiquettes des axes réduite
     cex.axis = 0.8)  # Définir les limites de l'axe Y de 1 à 3

# Ajouter la courbe pour x = b (correspondant à x_max) avec une ligne noire pointillée
lines(valeurs_y, f_xfixe_y[, 3], col = "black", lwd = 2, lty = 2)

# Colorier la région entre a (x_min) et b (x_max) en bleu transparent
polygon(c(valeurs_y, rev(valeurs_y)), 
        c(f_xfixe_y[, 1], rev(f_xfixe_y[, 3])), 
        col = rgb(0.1, 0.1, 0.8, 0.2), border = NA)  # Remplissage en bleu transparent

lines(valeurs_y, f_xfixe_y[, 2], col = "red", lwd = 2) # Ajouter la courbe pour la valeur médiane de x avec une ligne rouge
legend("bottomright", legend = c(paste("x =", a),  paste("x =", (a + b) / 2), paste("x =", b)),  col = c("black", "red", "black"), lwd = 2, lty = c(2, 1, 2)) # Ajouter une légende pour identifier les courbes

# Définir p2, la densité de probabilité uniforme pour Y
p2 <- function(y) 1  # Uniforme, donc constante entre c et d

# Définir p, la fonction de densité conjointe
p = function(x, y) p1(x) * p2(y)  # Densité conjointe p(x, y) = p1(x) * p2(y)

# Générer des échantillons pour X et Y pour la méthode de variance réduite
x_ech <- matrix(F1_inverse(runif(n_simulations * n_iterations)), ncol = n_iterations)  # Générer une matrice d'échantillons pour X selon p1 via F1_inverse
y_ech <- matrix(runif(n_simulations * n_iterations, c, d), ncol = n_iterations)  # Générer une matrice d'échantillons pour Y selon une loi uniforme entre c et d

# Calcul de f(x, y) et p(x, y) pour chaque échantillon tiré
valeurs_z_VR <- f_xy(x_ech, y_ech)  # Valeurs de f(x, y) pour chaque échantillon
valeurs_p_VR <- p(x_ech, y_ech)     # Valeurs de p(x, y) pour chaque échantillon

# Calculer l'évolution de la variance réduite (VR_evol)
VR_evol = cumsum(valeurs_z_VR[, 1] / valeurs_p_VR[, 1]) / (1:n_simulations)  # Moyenne cumulative pondérée par la fonction de densité p(x, y)

# Calcul de l'estimation de la variance réduite (VR) pour toutes les simulations
calculer_estimation_VR <- function(valeurs_z, valeurs_p) {
  mean(valeurs_z / valeurs_p)  # Moyenne des ratios f(x, y) / p(x, y)
}

# Appliquer la fonction de calcul sur chaque colonne de la matrice
VR_vecteur <- sapply(1:n_iterations, function(i) {
  calculer_estimation_VR(valeurs_z_VR[, i], valeurs_p_VR[, i])
})

#===============================================================================
# Méthode 4 : Méthode avec copules, marges méthode 2
#===============================================================================

# Fonction à intégrer avec rotation de 180 degrés
# Cette fonction applique une transformation à deux variables x et y
# en utilisant une fonction exponentielle et la fonction sinus,
# afin de créer une surface qui dépend des paramètres a, b, c et d.
f_xy_rotated <- function(x, y) {
  exp(sin((a + b - x) * (c + d - y)))  # Transformation des coordonnées
}

# Calcul des valeurs Z en utilisant 'outer' pour appliquer f_xy_rotated
# 'outer' crée une matrice de résultats en appliquant f_xy_rotated à chaque combinaison
# des valeurs x et y. Cela produit une surface Z correspondante.
valeurs_z_rotated <- outer(valeurs_x, valeurs_y, f_xy_rotated)

# Affichage du graphique 3D avec la fonction transformée
# On utilise 'persp3d' pour tracer la surface 3D avec les valeurs générées
persp3d(valeurs_x, valeurs_y, valeurs_z_rotated,
        col = "lightblue",                # Couleur de la surface
        xlab = "X",                       # Étiquette de l'axe X
        ylab = "Y",                       # Étiquette de l'axe Y
        zlab = "f(x, y)",                 # Étiquette de l'axe Z
        phi = 30,                         # Angle de la vue en élévation
        theta = 30,                       # Angle de la vue en rotation
        lit = TRUE,                       # Éclairage de la surface
        axes = TRUE)                      # Afficher les axes

# Paramètres et corrélation cible
tau_cible <- 0.5 # Valeur à ajuster pour la dépendance souhaitée entre X et Y

# Fonction pour minimiser la différence entre tau cible et tau généré
# Cette fonction calcule la différence absolue entre le tau calculé à partir de la copule
# et le tau_cible, afin de trouver le paramètre de dépendance optimal.
min_diff_tau <- function(theta) {
  copule <- claytonCopula(param = theta, dim = 2)  # Créer la copule de Clayton
  tau_calc <- tau(copule)  # Calcul du tau pour la copule de Clayton
  return(abs(tau_calc - tau_cible))  # Retourner la différence
}

# Optimisation pour trouver le theta optimal
# On utilise 'optimize' pour minimiser la différence entre tau calculé et tau cible
res <- optimize(min_diff_tau, interval = c(0.01, 10))  # Plage pour theta
theta_clayton <- res$minimum  # Extraire la valeur de theta optimale

# Affichage du theta optimal trouvé
cat("Theta optimal pour atteindre la dépendance cible : ", theta_clayton, "\n")

# Génération d'échantillons
n <- n_simulations * n_iterations  # Nombre total d'échantillons à générer

# Générer deux séries de variables uniformes indépendantes
U <- runif(n)  # Première variable uniforme
W <- runif(n)  # Variable auxiliaire pour générer la dépendance

# Calcul de V en utilisant la dépendance de la copule de Clayton
V <- (U^(-theta_clayton) * (W^(-theta_clayton / (1 + theta_clayton)) - 1) + 1)^(-1 / theta_clayton)

# Assemblage des paires (U, V)
echantillons_uv <- cbind(U, V)

# Transformation des échantillons (u, v) vers (x, y)
# La transformation convertit les marges uniformes [0, 1] en [a, b] pour X et [c, d] pour Y
x_ech <- matrix(qunif(echantillons_uv[, 1], min = a, max = b), ncol = n_iterations)  # Transforme U en X
y_ech <- matrix(qunif(echantillons_uv[, 2], min = c, max = d), ncol = n_iterations)  # Transforme V en Y

# Calcul de f_rotated(x, y) pour chaque échantillon
valeurs_z_copules <- f_xy_rotated(x_ech, y_ech)

# Calcul des probabilités marginales
u <- punif(x_ech, min = a, max = b)  # Probabilité pour X uniforme sur [a, b]
v <- punif(y_ech, min = c, max = d)  # Probabilité pour Y uniforme sur [c, d]

# Calcul des densités marginales
f_X <- dunif(x_ech, min = a, max = b)  # Densité de X
f_Y <- dunif(y_ech, min = c, max = d)  # Densité de Y

# Fonction de la densité de la copule de Clayton
c_clayton_density <- function(u, v, theta) {
  (1 + theta) * (u * v)^(-(1 + theta)) * (u^(-theta) + v^(-theta) - 1)^(-(2 + 1/theta))
}

# Calcul de la densité des échantillons selon Sklar
valeurs_p_copules <- c_clayton_density(u, v, theta_clayton) * f_X * f_Y

# Approximation de l'intégrale avec copules
# Cn_evol utilise la méthode des sommes cumulées pour approximer une intégrale
Cn_evol <- cumsum(valeurs_z_copules[, 1] / valeurs_p_copules[, 1]) / (1:n_simulations)

# Calcul de l'estimation de Cn pour toutes les simulations
# Fonction pour calculer l'estimation de Cn en utilisant les valeurs de Z et de P
calculer_estimation_copules <- function(valeurs_z, valeurs_p) {
  mean(valeurs_z / valeurs_p)  # Calcul de l'estimation
}

# Appliquer la fonction de calcul sur chaque colonne de la matrice
Cn_vecteur <- sapply(1:n_iterations, function(i) {
  calculer_estimation_copules(valeurs_z_copules[, i], valeurs_p_copules[, i])  # Calcul pour chaque itération
})

#===============================================================================
# Méthode 5 : Méthode avec copules, marges méthode 3
#===============================================================================

# Calculer les valeurs de f_xy_rotated(x, y_fixe) pour chaque combinaison de x et y
f_x_yfixe_rotated <- sapply(y_fixe, function(y) sapply(valeurs_x, function(x) f_xy_rotated(x, y)))  
# Création d'une matrice où chaque colonne correspond à f_xy_rotated(x, y) pour une valeur de y fixée

# Tracer la courbe pour y = c (correspondant à y_min) avec une ligne noire pointillée
plot(valeurs_x, f_x_yfixe_rotated[, 1], type = "l", col = "black", lwd = 2, lty = 2,
     xlab = "X", ylab = "f_rotated(X, Y)", 
     main = "Graphique de f_rotated(x, y) pour y : c, médiane, d",  # Titre plus explicite
     ylim = c(1, 3),
     cex.main = 0.8,   # Taille du titre réduite
     cex.lab = 0.8,    # Taille des étiquettes des axes réduite
     cex.axis = 0.8)   # Taille des étiquettes des axes réduite

# Ajouter la courbe pour y = d (correspondant à y_max) avec une ligne noire pointillée
lines(valeurs_x, f_x_yfixe_rotated[, 3], col = "black", lwd = 2, lty = 2)

# Colorier la région entre y = c (y_min) et y = d (y_max) en bleu transparent
polygon(c(valeurs_x, rev(valeurs_x)), 
        c(f_x_yfixe_rotated[, 1], rev(f_x_yfixe_rotated[, 3])), 
        col = rgb(0.1, 0.1, 0.8, 0.2), border = NA)  # Remplissage en bleu transparent

# Ajouter la courbe pour la valeur médiane de y avec une ligne rouge
lines(valeurs_x, f_x_yfixe_rotated[, 2], col = "red", lwd = 2)

# Calculer et tracer la fonction p1 en bleu
p1_rotated_valeurs <- 8 / 7 * (2.5 - 3 * valeurs_x)  # Définition de p1 en fonction de x
lines(valeurs_x, p1_rotated_valeurs, col = "blue", lwd = 2)  # Tracer la courbe de p1 en bleu

# Ajouter une légende pour identifier les courbes
legend("topright", 
       legend = c(paste("y =", c),          # Courbe pour y = c (y_min)
                  paste("y =", (c + d) / 2), # Courbe pour la valeur médiane de y
                  paste("y =", d),           # Courbe pour y = d (y_max)
                  "p1_rotated = 8/7 * (2.5 - 3x)"),  # Formule pour p1
       col = c("black", "red", "black", "blue"),  # Couleurs des courbes dans la légende
       lwd = 2, lty = c(2, 1, 2, 1))  # Styles de lignes dans la légende

# Afficher une grille sur le graphique
grid()

# Définition de p1_rotated, la densité de probabilité pour X
p1_rotated <- function(x) 8 / 7 * (2.5 - 3 * x)

# Fonction de répartition cumulée de p1_rotated
F1_rotated <- function(t) {
  return(4 / 7 * (5 * t - 3 * t^2))  # F1_rotated(t) correspond à la fonction de répartition cumulée de p1_rotated
}

# Inverse de la fonction de répartition F1_rotated (utile pour générer des échantillons de X selon p1_rotated)
F1_rotated_inverse <- function(z) {
  return((20 - sqrt(400 - 336 * z)) / 24)  # Résolution analytique de l'inverse de F1_rotated
}

# Histogramme des valeurs générées par la fonction inverse de F1_rotated avec une distribution uniforme
hist(F1_rotated_inverse(runif(10000)), 
     col = rgb(0.1, 0.5, 0.8, 0.7),  # Couleur des barres avec transparence
     border = "white",  # Couleur de la bordure des barres
     xlab = "Valeurs de X",  # Étiquette de l'axe X
     ylab = "Fréquence",  # Étiquette de l'axe Y
     main = paste("Histogramme des échantillons de X", "\nselon la densité de proba 1"),  # Titre avec saut de ligne
     xlim = c(0, 0.5),  # Limites de l'axe X (ajuster selon votre intervalle)
     ylim = c(0, 1400),  # Limites de l'axe Y (ajuster selon votre distribution)
     las = 1,  # Orientation des étiquettes de l'axe Y
     freq = TRUE,
     cex.main = 0.8,   # Taille du titre réduite
     cex.lab = 0.8,    # Taille des étiquettes des axes réduite
     cex.axis = 0.8)  # Fréquence ou densité (TRUE pour fréquence)

# Calculer les valeurs de f_xy_rotated(x, y_fixe) pour chaque combinaison de x et y
f_xfixe_y_rotated <- sapply(x_fixe, function(x) sapply(valeurs_y, function(y) f_xy_rotated(x, y)))  
# Création d'une matrice où chaque colonne correspond à f_xy_rotated(x, y) pour une valeur de x fixée

# Tracer la courbe pour x = a (correspondant à x_min) avec une ligne noire pointillée
plot(valeurs_y, f_xfixe_y_rotated[, 1], type = "l", col = "black", lwd = 2, lty = 2,
     xlab = "Y", ylab = "f_rotated(X, Y)", 
     main = "Graphique de f_rotated(x, y) pour x : a, médiane, b",  # Titre plus clair
     ylim = c(1, 3),
     cex.main = 0.8,   # Taille du titre réduite
     cex.lab = 0.8,    # Taille des étiquettes des axes réduite
     cex.axis = 0.8)  # Définir les limites de l'axe Y de 1 à 3

# Ajouter la courbe pour x = b (correspondant à x_max) avec une ligne noire pointillée
lines(valeurs_y, f_xfixe_y_rotated[, 3], col = "black", lwd = 2, lty = 2)

# Colorier la région entre a (x_min) et b (x_max) en bleu transparent
polygon(c(valeurs_y, rev(valeurs_y)), 
        c(f_xfixe_y_rotated[, 1], rev(f_xfixe_y_rotated[, 3])), 
        col = rgb(0.1, 0.1, 0.8, 0.2), border = NA)  # Remplissage en bleu transparent

lines(valeurs_y, f_xfixe_y_rotated[, 2], col = "red", lwd = 2) # Ajouter la courbe pour la valeur médiane de x avec une ligne rouge
legend("bottomright", legend = c(paste("x =", a),  paste("x =", (a + b) / 2), paste("x =", b)),  col = c("black", "red", "black"), lwd = 2, lty = c(2, 1, 2)) # Ajouter une légende pour identifier les courbes

# Définir p2_rotated, la densité de probabilité uniforme pour Y
p2_rotated <- function(y) 1  # Uniforme, donc constante entre c et d

# Définir p_rotated, la fonction de densité conjointe
p_rotated = function(x, y) p1_rotated(x) * p2_rotated(y)  # Densité conjointe p_rotated(x, y) = p1_rotated(x) * p2_rotated(y)

# Générer deux séries de variables uniformes indépendantes
U <- runif(n)  # Première variable uniforme
W <- runif(n)  # Variable auxiliaire pour générer la dépendance

# Calcul de V en utilisant la dépendance de la copule de Clayton
V <- (U^(-theta_clayton) * (W^(-theta_clayton / (1 + theta_clayton)) - 1) + 1)^(-1 / theta_clayton)

# Assemblage des paires (U, V)
echantillons_uv <- cbind(U, V)

# Transformation des échantillons (u, v) vers (x, y)
# La transformation convertit les marges uniformes [0, 1] en [a, b] pour X et [c, d] pour Y
x_ech <- matrix(F1_rotated_inverse(echantillons_uv[, 1]), ncol = n_iterations)  # Transforme U en X
y_ech <- matrix(qunif(echantillons_uv[, 2], min = c, max = d), ncol = n_iterations)  # Transforme V en Y

# Calcul de f_rotated(x, y) pour chaque échantillon
valeurs_z_copules2 <- f_xy_rotated(x_ech, y_ech)

# Calcul des probabilités marginales
u <- F1_rotated(x_ech)  # Probabilité pour X uniforme sur [a, b]
v <- punif(y_ech, min = c, max = d)  # Probabilité pour Y uniforme sur [c, d]

# Calcul des densités marginales
f_X <- p1_rotated(x_ech)  # Densité de X
f_Y <- dunif(y_ech, min = c, max = d)  # Densité de Y

# Calcul de la densité des échantillons selon Sklar
valeurs_p_copules2 <- c_clayton_density(u, v, theta_clayton) * f_X * f_Y

# Approximation de l'intégrale avec copules
# Cn_evol utilise la méthode des sommes cumulées pour approximer une intégrale
Cn_2_evol <- cumsum(valeurs_z_copules2[, 1] / valeurs_p_copules2[, 1]) / (1:n_simulations)

# Appliquer la fonction de calcul sur chaque colonne de la matrice
Cn_vecteur2 <- sapply(1:n_iterations, function(i) {
  calculer_estimation_copules(valeurs_z_copules2[, i], valeurs_p_copules2[, i])  # Calcul pour chaque itération
})

#===============================================================================
# Méthode 6 : Méthode Quasi-Monte Carlo (avec séquences de Halton)
#===============================================================================

# Génération de points de Halton pour le nombre total de simulations et d'itérations
halton_points <- halton(n_simulations * n_iterations, dim = 2)

# Mise à l'échelle des points de Halton pour les limites définies
# x_ech et y_ech sont des matrices où chaque colonne correspond à une itération
x_ech <- matrix(a + (b - a) * halton_points[, 1], ncol = n_iterations)
y_ech <- matrix(c + (d - c) * halton_points[, 2], ncol = n_iterations)

# Calcul de f(x, y) pour chaque échantillon tiré
valeurs_z_halton <- f_xy(x_ech, y_ech)

# Approximation de l'intégrale par Quasi-Monte Carlo
# aire_rectangle est le volume du domaine d'intégration
aire_rectangle <- (b - a) * (d - c)
Hn_evol <- cumsum(valeurs_z_halton[, 1]) / (1:n_simulations) * aire_rectangle

# Fonction pour calculer l'estimation de l'intégrale
calculer_estimation_qmc <- function(valeurs) {
  mean(valeurs) * aire_rectangle  # Retourne la moyenne des valeurs multipliée par le volume
}

# Calcul des estimations pour toutes les simulations
Hn_vecteur <- apply(valeurs_z_halton, MARGIN = 2, FUN = calculer_estimation_qmc)

#===============================================================================
# Calcul des risques quadratiques
#===============================================================================

# Estimation de l'intégrale par un calcul numérique
# Utilisation de la fonction integral2 pour calculer la valeur exacte de l'intégrale de la fonction f_xy
I <- integral2(f_xy, a, b, c, d)$Q

# Calcul des risques quadratiques pour chaque méthode
# RQ_Tn : Risque quadratique pour la méthode Volume (Tn)
# RQ_Sn : Risque quadratique pour la méthode Monte Carlo basique (Sn)
# RQ_VR : Risque quadratique pour la méthode de Variance Réduite (VR)
# RQ_Cn : Risque quadratique pour la méthode Monte Carlo avec Copules (Cn)
# RQ_Cn2 : Risque quadratique pour la méthode Monte Carlo avec Copules (Cn2)
# RQ_Hn : Risque quadratique pour la méthode Quasi-Monte Carlo (Hn)
RQ_Tn <- (mean(Tn_vecteur) - I)^2 + sd(Tn_vecteur)^2
RQ_Sn <- (mean(Sn_vecteur) - I)^2 + sd(Sn_vecteur)^2
RQ_VR <- (mean(VR_vecteur) - I)^2 + sd(VR_vecteur)^2
RQ_Cn <- (mean(Cn_vecteur) - I)^2 + sd(Cn_vecteur)^2
RQ_Cn2 <- (mean(Cn_vecteur2) - I)^2 + sd(Cn_vecteur2)^2 
RQ_Hn <- (mean(Hn_vecteur) - I)^2 + sd(Hn_vecteur)^2

#===============================================================================
# Résultats et visualisation
#===============================================================================

# Affichage des résultats des estimations d'intégrale pour chaque méthode
cat("Valeur approximée par la méthode Volume (Tn): ", mean(Tn_vecteur), "\n")
cat("Valeur approximée par la méthode Monte Carlo basique (Sn): ", mean(Sn_vecteur), "\n")
cat("Valeur approximée par la méthode Variance Réduite (VR): ", mean(VR_vecteur), "\n")
cat("Valeur approximée par la méthode Monte Carlo avec Copules (Cn): ", mean(Cn_vecteur), "\n")
cat("Valeur approximée par la méthode Monte Carlo avec Copules (Cn2): ", mean(Cn_vecteur2), "\n")
cat("Valeur approximée par la méthode Quasi-Monte Carlo (Hn): ", mean(Hn_vecteur), "\n")

# Affichage des risques quadratiques pour chaque méthode
cat("Risque quadratique pour Tn: ", RQ_Tn, "\n")
cat("Risque quadratique pour Sn: ", RQ_Sn, "\n")
cat("Risque quadratique pour VR: ", RQ_VR, "\n")
cat("Risque quadratique pour Cn: ", RQ_Cn, "\n")
cat("Risque quadratique pour Cn2: ", RQ_Cn2, "\n")
cat("Risque quadratique pour Hn: ", RQ_Hn, "\n")

# Graphique de comparaison des méthodes
# Création d'un graphique avec les évolutions des valeurs estimées pour chaque méthode
plot(1:n_simulations, Tn_evol, type = "l", col = "blue", lwd = 2,
     xlab = "Nombre de simulations", ylab = "Valeur estimée de l'intégrale",
     main = "Comparaison des méthodes d'estimation de l'intégrale",
     xlim = c(0, 1000),
     ylim = range(c(Tn_evol, Sn_evol, VR_evol, Cn_evol, Cn_2_evol, Hn_evol, I)),
     las = 1,
     cex.main = 0.8,
     cex.lab = 0.8,
     cex.axis = 0.8)

# Ajout des lignes pour chaque méthode
lines(1:n_simulations, Sn_evol, col = "purple", lwd = 2)
lines(1:n_simulations, VR_evol, col = "green", lwd = 2)
lines(1:n_simulations, Cn_evol, col = "pink", lwd = 2)
lines(1:n_simulations, Cn_2_evol, col = "orange", lwd = 2)
lines(1:n_simulations, Hn_evol, col = "cyan", lwd = 2)

# Ajout d'une ligne horizontale représentant la valeur exacte de l'intégrale
abline(h = I, col = "red", lty = 2, lwd = 2)

# Création de la légende pour le graphique
legend("topright", legend = c("Méthode Volume (Tn)", "Monte Carlo Basique (Sn)", 
                              "Variance Réduite (VR)", "Méthode des copules (Cn)", 
                              "Méthode des copules (Cn2)", "Méthode Quasi-Monte Carlo (Hn)", 
                              "Valeur exacte"),
       col = c("blue", "purple", "green", "pink", "orange", "cyan", "red"),
       lty = c(1, 1, 1, 1, 1, 1, 2), lwd = 2, bty = "n", cex = 0.8)
grid()

# Fonction pour afficher l'histogramme d'un vecteur de valeurs avec des limites x spécifiques
# Affiche l'histogramme des estimations pour visualiser la distribution des valeurs
afficher_histogramme <- function(vecteur, titre, I) {
  xlim_individuel <- range(c(vecteur, I))  # Limites x spécifiques au vecteur
  hist(vecteur, main = paste("Distribution des estimations", "\n(", titre, ")"), 
       xlab = "Estimation", 
       col = "lightblue", 
       border = "darkblue", 
       breaks = 30,  
       xlim = xlim_individuel, 
       las = 1,  
       cex.main = 0.8)
  abline(v = I, col = "red", lty = 2, lwd = 2)  # Ligne verticale pour la valeur de I
  grid()
}

# Organisation des histogrammes
par(mfrow = c(3, 2))  # Organisation en 3x2 pour les graphiques

# Affichage des histogrammes avec des limites x spécifiques pour chaque méthode
afficher_histogramme(Tn_vecteur, "Méthode Volume", I)
afficher_histogramme(Sn_vecteur, "Méthode MC Basique", I)
afficher_histogramme(VR_vecteur, "Méthode Variance Réduite", I)
afficher_histogramme(Cn_vecteur, "Méthode Monte Carlo avec Copules (Cn)", I)
afficher_histogramme(Cn_vecteur2, "Méthode Monte Carlo avec Copules (Cn2)", I)
afficher_histogramme(Hn_vecteur, "Méthode Quasi-Monte Carlo (Hn)", I) 

# Réinitialisation de l'affichage à une seule fenêtre graphique
par(mfrow = c(1, 1))  # Réinitialisation de l'affichage
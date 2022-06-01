# ---------------------------------
# Knockoffs vs Knockoffs+
# ---------------------------------
source("Functions.R")

# Number of selected variables
# A = 1
a1 <- read.table('Simulation Study 2/selbin_a1_q1.txt')
a2 <- read.table('Simulation Study 2/selbin_a1_q2.txt')
a3 <- read.table('Simulation Study 2/selbin_a1_q3.txt')
a4 <- read.table('Simulation Study 2/selbin_a1_q4.txt')
a5 <- read.table('Simulation Study 2/selbin_a1_q5.txt')
plusa1 <- read.table('Simulation Study 2/selbin_pa1_q1.txt')
plusa2 <- read.table('Simulation Study 2/selbin_pa1_q2.txt')
plusa3 <- read.table('Simulation Study 2/selbin_pa1_q3.txt')
plusa4 <- read.table('Simulation Study 2/selbin_pa1_q4.txt')
plusa5 <- read.table('Simulation Study 2/selbin_pa1_q5.txt')

a1sel <- comp_nr_selected(a1, 1)
a2sel <- comp_nr_selected(a2, 1)
a3sel <- comp_nr_selected(a3, 1)
a4sel <- comp_nr_selected(a4, 1)
a5sel <- comp_nr_selected(a5, 1)
plusa1sel <- comp_nr_selected(plusa1, 1)
plusa2sel <- comp_nr_selected(plusa2, 1)
plusa3sel <- comp_nr_selected(plusa3, 1)
plusa4sel <- comp_nr_selected(plusa4, 1)
plusa5sel <- comp_nr_selected(plusa5, 1)

a1_mean <- comp_mean_n(c(a1sel, a2sel, a3sel, a4sel, a5sel, 
                         plusa1sel, plusa2sel, plusa3sel, plusa4sel, plusa5sel), 20)




# A = 5
n1 <- read.table('Simulation Study 2/selbin_a5_q1.txt')
n2 <- read.table('Simulation Study 2/selbin_a5_q2.txt')
n3 <- read.table('Simulation Study 2/selbin_a5_q3.txt')
n4 <- read.table('Simulation Study 2/selbin_a5_q4.txt')
n5 <- read.table('Simulation Study 2/selbin_a5_q5.txt')
plusn1 <- read.table('Simulation Study 2/selbin_pa5_q1.txt')
plusn2 <- read.table('Simulation Study 2/selbin_pa5_q2.txt')
plusn3 <- read.table('Simulation Study 2/selbin_pa5_q3.txt')
plusn4 <- read.table('Simulation Study 2/selbin_pa5_q4.txt')
plusn5 <- read.table('Simulation Study 2/selbin_pa5_q5.txt')


n1sel <- comp_nr_selected(n1, 1)
n2sel <- comp_nr_selected(n2, 1)
n3sel <- comp_nr_selected(n3, 1)
n4sel <- comp_nr_selected(n4, 1)
n5sel <- comp_nr_selected(n5, 1)
plusn1sel <- comp_nr_selected(plusn1, 1)
plusn2sel <- comp_nr_selected(plusn2, 1)
plusn3sel <- comp_nr_selected(plusn3, 1)
plusn4sel <- comp_nr_selected(plusn4, 1)
plusn5sel <- comp_nr_selected(plusn5, 1)


a5_mean <- comp_mean_n(c(n1sel, n2sel, n3sel, n4sel, n5sel, 
                         plusn1sel, plusn2sel, plusn3sel, plusn4sel, plusn5sel), 20)




# A = 10
w1 <- read.table('Simulation Study 2/selbin_a10_q1.txt')
w2 <- read.table('Simulation Study 2/selbin_a10_q2.txt')
w3 <- read.table('Simulation Study 2/selbin_a10_q3.txt')
w4 <- read.table('Simulation Study 2/selbin_a10_q4.txt')
w5 <- read.table('Simulation Study 2/selbin_a10_q5.txt')
plusw1 <- read.table('Simulation Study 2/selbin_pa10_q1.txt')
plusw2 <- read.table('Simulation Study 2/selbin_pa10_q2.txt')
plusw3 <- read.table('Simulation Study 2/selbin_pa10_q3.txt')
plusw4 <- read.table('Simulation Study 2/selbin_pa10_q4.txt')
plusw5 <- read.table('Simulation Study 2/selbin_pa10_q5.txt')


w1sel <- comp_nr_selected(w1, 1)
w2sel <- comp_nr_selected(w2, 1)
w3sel <- comp_nr_selected(w3, 1)
w4sel <- comp_nr_selected(w4, 1)
w5sel <- comp_nr_selected(w5, 1)
plusw1sel <- comp_nr_selected(plusw1, 1)
plusw2sel <- comp_nr_selected(plusw2, 1)
plusw3sel <- comp_nr_selected(plusw3, 1)
plusw4sel <- comp_nr_selected(plusw4, 1)
plusw5sel <- comp_nr_selected(plusw5, 1)


a10_mean <- comp_mean_n(c(w1sel, w2sel, w3sel, w4sel, w5sel, 
                         plusw1sel, plusw2sel, plusw3sel, plusw4sel, plusw5sel), 20)


# A = 20
z1 <- read.table('Simulation Study 2/selbin_a20_q1.txt')
z2 <- read.table('Simulation Study 2/selbin_a20_q2.txt')
z3 <- read.table('Simulation Study 2/selbin_a20_q3.txt')
z4 <- read.table('Simulation Study 2/selbin_a20_q4.txt')
z5 <- read.table('Simulation Study 2/selbin_a20_q5.txt')
plusz1 <- read.table('Simulation Study 2/selbin_pa20_q1.txt')
plusz2 <- read.table('Simulation Study 2/selbin_pa20_q2.txt')
plusz3 <- read.table('Simulation Study 2/selbin_pa20_q3.txt')
plusz4 <- read.table('Simulation Study 2/selbin_pa20_q4.txt')
plusz5 <- read.table('Simulation Study 2/selbin_pa20_q5.txt')


z1sel <- comp_nr_selected(z1, 1)
z2sel <- comp_nr_selected(z2, 1)
z3sel <- comp_nr_selected(z3, 1)
z4sel <- comp_nr_selected(z4, 1)
z5sel <- comp_nr_selected(z5, 1)
plusz1sel <- comp_nr_selected(plusz1, 1)
plusz2sel <- comp_nr_selected(plusz2, 1)
plusz3sel <- comp_nr_selected(plusz3, 1)
plusz4sel <- comp_nr_selected(plusz4, 1)
plusz5sel <- comp_nr_selected(plusz5, 1)

a20_mean <- comp_mean_n(c(z1sel, z2sel, z3sel, z4sel, z5sel, 
                          plusz1sel, plusz2sel, plusz3sel, plusz4sel, plusz5sel), 20)



# A = 30
x1 <- read.table('Simulation Study 2/selbin_a30_q1.txt')
x2 <- read.table('Simulation Study 2/selbin_a30_q2.txt')
x3 <- read.table('Simulation Study 2/selbin_a30_q3.txt')
x4 <- read.table('Simulation Study 2/selbin_a30_q4.txt')
x5 <- read.table('Simulation Study 2/selbin_a30_q5.txt')
plusx1 <- read.table('Simulation Study 2/selbin_pa30_q1.txt')
plusx2 <- read.table('Simulation Study 2/selbin_pa30_q2.txt')
plusx3 <- read.table('Simulation Study 2/selbin_pa30_q3.txt')
plusx4 <- read.table('Simulation Study 2/selbin_pa30_q4.txt')
plusx5 <- read.table('Simulation Study 2/selbin_pa30_q5.txt')


x1sel <- comp_nr_selected(x1, 1)
x2sel <- comp_nr_selected(x2, 1)
x3sel <- comp_nr_selected(x3, 1)
x4sel <- comp_nr_selected(x4, 1)
x5sel <- comp_nr_selected(x5, 1)
plusx1sel <- comp_nr_selected(plusx1, 1)
plusx2sel <- comp_nr_selected(plusx2, 1)
plusx3sel <- comp_nr_selected(plusx3, 1)
plusx4sel <- comp_nr_selected(plusx4, 1)
plusx5sel <- comp_nr_selected(plusx5, 1)


a30_mean <- comp_mean_n(c(x1sel, x2sel, x3sel, x4sel, x5sel, 
                          plusx1sel, plusx2sel, plusx3sel, plusx4sel, plusx5sel), 20)



# A = 40
y1 <- read.table('Simulation Study 2/selbin_a40_q1.txt')
y2 <- read.table('Simulation Study 2/selbin_a40_q2.txt')
y3 <- read.table('Simulation Study 2/selbin_a40_q3.txt')
y4 <- read.table('Simulation Study 2/selbin_a40_q4.txt')
y5 <- read.table('Simulation Study 2/selbin_a40_q5.txt')
plusy1 <- read.table('Simulation Study 2/selbin_pa40_q1.txt')
plusy2 <- read.table('Simulation Study 2/selbin_pa40_q2.txt')
plusy3 <- read.table('Simulation Study 2/selbin_pa40_q3.txt')
plusy4 <- read.table('Simulation Study 2/selbin_pa40_q4.txt')
plusy5 <- read.table('Simulation Study 2/selbin_pa40_q5.txt')


y1sel <- comp_nr_selected(y1, 1)
y2sel <- comp_nr_selected(y2, 1)
y3sel <- comp_nr_selected(y3, 1)
y4sel <- comp_nr_selected(y4, 1)
y5sel <- comp_nr_selected(y5, 1)
plusy1sel <- comp_nr_selected(plusy1, 1)
plusy2sel <- comp_nr_selected(plusy2, 1)
plusy3sel <- comp_nr_selected(plusy3, 1)
plusy4sel <- comp_nr_selected(plusy4, 1)
plusy5sel <- comp_nr_selected(plusy5, 1)


a40_mean <- comp_mean_n(c(y1sel, y2sel, y3sel, y4sel, y5sel, 
                          plusy1sel, plusy2sel, plusy3sel, plusy4sel, plusy5sel), 20)



# A = 50
h1 <- read.table('Simulation Study 2/selbin_a40_q1.txt')
h2 <- read.table('Simulation Study 2/selbin_a40_q2.txt')
h3 <- read.table('Simulation Study 2/selbin_a40_q3.txt')
h4 <- read.table('Simulation Study 2/selbin_a40_q4.txt')
h5 <- read.table('Simulation Study 2/selbin_a40_q5.txt')
plush1 <- read.table('Simulation Study 2/selbin_pa40_q1.txt')
plush2 <- read.table('Simulation Study 2/selbin_pa40_q2.txt')
plush3 <- read.table('Simulation Study 2/selbin_pa40_q3.txt')
plush4 <- read.table('Simulation Study 2/selbin_pa40_q4.txt')
plush5 <- read.table('Simulation Study 2/selbin_pa40_q5.txt')


h1sel <- comp_nr_selected(h1, 1)
h2sel <- comp_nr_selected(h2, 1)
h3sel <- comp_nr_selected(h3, 1)
h4sel <- comp_nr_selected(h4, 1)
h5sel <- comp_nr_selected(h5, 1)
plush1sel <- comp_nr_selected(plush1, 1)
plush2sel <- comp_nr_selected(plush2, 1)
plush3sel <- comp_nr_selected(plush3, 1)
plush4sel <- comp_nr_selected(plush4, 1)
plush5sel <- comp_nr_selected(plush5, 1)

a50_mean <- comp_mean_n(c(h1sel, h2sel, h3sel, h4sel, h5sel, 
                          plush1sel, plush2sel, plush3sel, plush4sel, plush5sel), 20)



# ---------------------------
# Data Frames for plotting
# ---------------------------

df <- data.frame(nsel = c(a1_mean, a5_mean, a10_mean, a20_mean, a30_mean, a40_mean, a50_mean),
                 q = as.factor(c(0.05, 0.1, 0.2, 0.3, 0.5)),
                 amp = c(rep(1, 10), rep(5, 10), rep(10, 10), rep(20, 10), rep(30, 10), rep(40, 10), rep(50, 10)),
                 Method = c(rep("Knockoffs", 5), rep("Knockoffs+", 5)))




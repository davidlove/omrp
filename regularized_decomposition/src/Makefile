# Makefile for rd

CC=g++
CFLAGS=-ggdb -O3 -Wall
# CFLAGS=-O3 -Wall
LDFLAGS=-lm
SOURCES=changelp.cpp cl_list.cpp config.cpp crash.cpp cust_opt.cpp \
 	 determlp.cpp error.cpp gnu_fix1.cpp gnu_fix2.cpp gnu_fix3.cpp \
	 gnu_fix4.cpp history.cpp init_bas.cpp invaux.cpp inverse.cpp \
	 invfact.cpp invsolve.cpp invupd.cpp io_dit.cpp lexer.cpp lp_sol.cpp \
	 main.cpp main_arg.cpp makelab.cpp mast_sol.cpp memblock.cpp \
	 mps_lp.cpp num_elim.cpp option.cpp parsemps.cpp parsespc.cpp \
	 parsstoc.cpp periods.cpp postsolv.cpp pp_integ.cpp pp_primi.cpp \
	 pp_simpl.cpp pres_act.cpp presolve.cpp print.cpp q1.cpp q2.cpp \
	 q3.cpp q4.cpp q5.cpp q6.cpp rand01.cpp rd_solv.cpp rd_sublp.cpp \
	 read_act.cpp read_lp.cpp read_tim.cpp rm_infea.cpp row_anal.cpp \
	 sc_tree.cpp scenario.cpp scentree.cpp sing_col.cpp smplx_lp.cpp \
	 sol_lab.cpp solution.cpp solv_lp.cpp solver1.cpp solver2.cpp \
	 solver3.cpp solver4.cpp solvpiv.cpp solvpric.cpp solvtol.cpp \
	 sort_lab.cpp start_pt.cpp std_math.cpp stochsol.cpp strdup.cpp \
	 sub_man.cpp time_cnt.cpp vec_pool.cpp writ_mps.cpp
OBJS=	 $(SOURCES:.cpp=.o)

INVS=inverse.o invfact.o invupd.o
EXECUTABLE=../decomp

$(EXECUTABLE): $(OBJS)
	$(CC) $(LDFLAGS) $(OBJS) -o $@

${OBJS}: error.h memblock.h sptr_deb.h stdtype.h smartdcl.h smartptr.h

cl_list.o gnu_fix1.o io_dit.o: sptr_ndb.h myalloc.h

crash.o cust_opt.o determlp.o history.o init_bas.o $(INVS) io_dit.o lexer.o lp_sol.o main_arg.o main.o: smplx_lp.h solv_lp.h std_tmpl.h lp_codes.h mps_lp.h compile.h simplex.h sort_lab.h std_math.h presolve.h pre_code.h cl_list.h scenario.h option.h my_defs.h solvcode.h parsespc.h solution.h

crash.o main_arg.o main.o: print.h

gnu_fix1.o gnu_fix3.o gnu_fix4.o $(INVS) io_dit.o main.o: work_vec.h vec_pool.h

$(INVS) main.o: inverse.h invaux.h

clean:
	rm -f *.o

# I'm up to makelab.o in creating this full makefile

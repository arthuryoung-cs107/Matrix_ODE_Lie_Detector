# file names

EXHIB:= LD_aux LD_ode LD_io LD_framework LD_function_space LD_framework
EXHIBGEN:= LD_ode_system LD_integrators

include config/includes.mak

# Lie detector object files

EXHIBOBJS:= $(addprefix $(LD_EXHIB_DIR), $(addsuffix .o, $(EXHIB)))
EXHIBGENOBJS:= $(addprefix $(LD_EXHIB_DIR), $(addsuffix .o, $(EXHIBGEN)))

# final executable object files
LDEXHIBOBJS:= $(EXHIBOBJS) #"Lie Detector reconstruction"
LDEXHIBGENOBJS:= $(EXHIBOBJS) $(EXHIBGENOBJS) #"Lie Detector reconstruction"

# Lie detector exhibition rules
$(LD_EXHIB_DIR)%.o: $(LD_EXHIB_SRC)%.cc | $(LD_EXHIB_DIR)
	$(CXX) $(IDIR_EXHIB) $(CFLAGS_EXHIB) -c $< -o $@

denoise_demo: $(TEST_SRC)denoising_demo.cc $(LDEXHIBOBJS)
	$(CXX) $(IDIR_EXHIB) $(CFLAGS_EXHIB) $(LINK) $^ $(LIBS_EXHIB) -o $@

gendata_demo: $(TEST_SRC)generate_data_demo.cc $(LDEXHIBGENOBJS)
	$(CXX) $(IDIR_EXHIB) $(CFLAGS_EXHIB) $(LINK) $^ $(LIBS_EXHIB) -o $@

$(LD_EXHIB_DIR):
	mkdir -p $@

clean_:
	rm -f *_demo

clean_LD_EXHIB:
	rm -f $(LD_EXHIB_DIR)*.o

clean: clean_

clean_all: clean_ clean_LD_EXHIB

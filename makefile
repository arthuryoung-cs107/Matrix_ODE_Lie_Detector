# file names

EXHIB:= LD_ode LD_function_space LD_parameter_space LD_framework
EXHIBDEC:= matrix_Lie_detector
EXHIBGEN:= LD_integrators
EXHIBCOK:= LD_cokernals

include config/includes.mak

# Lie detector object files
EXHIBOBJS:= $(addprefix $(LD_EXHIB_DIR), $(addsuffix .o, $(EXHIB)))
EXHIBDECOBJS:= $(addprefix $(LD_EXHIB_DIR), $(addsuffix .o, $(EXHIBDEC)))
EXHIBGENOBJS:= $(addprefix $(LD_EXHIB_DIR), $(addsuffix .o, $(EXHIBGEN)))
EXHIBCOKOBJS:= $(addprefix $(LD_EXHIB_DIR), $(addsuffix .o, $(EXHIBCOK)))

# final executable object files
LDEXHIBOBJS:= $(EXHIBOBJS)
LDEXHIBDECOBJS:= $(EXHIBOBJS) $(EXHIBDECOBJS)
LDEXHIBGENOBJS:= $(EXHIBOBJS) $(EXHIBGENOBJS)
LDEXHIBRECOBJS:= $(EXHIBOBJS) $(EXHIBGENOBJS) $(EXHIBDECOBJS)
LDEXHIBCOKOBJS:= $(EXHIBOBJS) $(EXHIBGENOBJS) $(EXHIBDECOBJS) $(EXHIBCOKOBJS)

# Lie detector exhibition rules
$(LD_EXHIB_DIR)%.o: $(LD_EXHIB_SRC)%.cc | $(LD_EXHIB_DIR)
	$(CXX) $(IDIR_EXHIB) $(CFLAGS_EXHIB) -c $< -o $@

encode_demo: $(TEST_SRC)encoding_demo.cc $(LDEXHIBOBJS)
	$(CXX) $(IDIR_EXHIB) $(CFLAGS_EXHIB) $(LINK) $^ $(LIBS_EXHIB) -o $@

decode_demo: $(TEST_SRC)decoding_demo.cc $(LDEXHIBDECOBJS)
	$(CXX) $(IDIR_EXHIB) $(CFLAGS_EXHIB) $(LINK) $^ $(LIBS_EXHIB) -o $@

recon_demo: $(TEST_SRC)reconstruction_demo.cc $(LDEXHIBRECOBJS)
	$(CXX) $(IDIR_EXHIB) $(CFLAGS_EXHIB) $(LINK) $^ $(LIBS_EXHIB) -o $@

gendata_demo: $(TEST_SRC)generate_data_demo.cc $(LDEXHIBGENOBJS)
	$(CXX) $(IDIR_EXHIB) $(CFLAGS_EXHIB) $(LINK) $^ $(LIBS_EXHIB) -o $@

full_demo: $(TEST_SRC)full_demo.cc $(LDEXHIBRECOBJS)
	$(CXX) $(IDIR_EXHIB) $(CFLAGS_EXHIB) $(LINK) $^ $(LIBS_EXHIB) -o $@

transfer_demo: $(TEST_SRC)transfer_learning_demo.cc $(LDEXHIBCOKOBJS)
	$(CXX) $(IDIR_EXHIB) $(CFLAGS_EXHIB) $(LINK) $^ $(LIBS_EXHIB) -o $@

cokernal_demo: $(TEST_SRC)cokernal_demo.cc $(LDEXHIBCOKOBJS)
	$(CXX) $(IDIR_EXHIB) $(CFLAGS_EXHIB) $(LINK) $^ $(LIBS_EXHIB) -o $@

denoise_demo: $(TEST_SRC)denoise_demo.cc $(LDEXHIBCOKOBJS)
	$(CXX) $(IDIR_EXHIB) $(CFLAGS_EXHIB) $(LINK) $^ $(LIBS_EXHIB) -o $@

$(LD_EXHIB_DIR):
	mkdir -p $@

clean_:
	rm -f *_demo
	rm -f *_scratch

clean_debugging:
	rm -f *.lddat

clean_LD_EXHIB:
	rm -f $(LD_EXHIB_DIR)*.o

clean_data_directory:
	rm -f ./data_directory/*.lddat

clean: clean_

clean_all: clean_ clean_LD_EXHIB

clean_ALL: clean_all clean_data_directory clean_debugging

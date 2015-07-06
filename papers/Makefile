# typing 'make' here will build papers in all subdirectories

PAPERS = DNS_Reionization_I \
         DNS_Reionization_II 


all : 
	for dir in $(PAPERS); do \
	  $(MAKE) -C $$dir; \
	done

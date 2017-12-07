# Alexander Drabent, December 2017
# make losoto v2.0 parset file

import os

def plugin_main(args, **kwargs):

	steps = kwargs['steps'].lstrip('[').rstrip(']').replace(' ','').split(',')
	
	parset_dict = {}
	parset_dict['global'] = {}
	for step in steps:
		parset_dict[step.strip()] = {}
		pass
    
	for key in kwargs.keys():
		keyword = key.split('.')[0]
		option  = key.split('.')[-1]
		value   = kwargs[key]
		if keyword in steps:
			parset_dict[keyword][option] = value
			pass
		elif keyword == 'global':
			parset_dict[keyword][option] = value
		pass
          
	parset_file = open(kwargs['filename'], 'w')
	for option in parset_dict['global']:
		parset_file.write(option + ' = ' + parset_dict['global'][option] + '\n')
		pass

	for step in steps:
		parset_file.write('[' + step + ']\n')
		if step in parset_dict.keys():
			for option in parset_dict[step]: 
				parset_file.write(option.replace(';','.') + ' = ' + parset_dict[step][option] + '\n')
				pass
			pass
		pass
	parset_file.close()
       
	return 0
	pass
  

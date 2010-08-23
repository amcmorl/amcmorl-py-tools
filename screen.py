from subprocess import Popen, PIPE

def get_screen_resolution():
    output = Popen(['xwininfo', '-root'], stdout=PIPE).communicate()[0]
    output = [x.strip() for x in output.split('\n')]
    as_dict = {}
    for line in output:
        if ':' in line:
            key, value = [x.strip() for x in line.split(':', 1)]
            as_dict[key] = value
    return int(as_dict['Width']), int(as_dict['Height'])
        
def get_screen_resolution2():
    output = Popen(['xwininfo', '-root'], stdout=PIPE).communicate()[0]
    output = [x.strip() for x in output.split('\n')]
    as_dict = {}
    for line in output:
        if 'Width:' in line:
            width = int(line.split(':', 1)[1])
        elif 'Height' in line:
            height = int(line.split(':', 1)[1])
    return width, height


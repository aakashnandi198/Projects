
import os

def getCleanName(projectName):
    return (" ").join([x for x in projectName.split('_')])


if __name__ == '__main__':
    projects = [x for x in os.scandir("./")if x.is_dir() and not x.name.startswith('.')]
    cnt = 1
    print('\n-------------Projects Available--------------')
    for project in os.scandir("./"):
        if project.is_dir() and not project.name.startswith('.'):
            print(str(cnt)+") "+getCleanName(project.name))
            cnt += 1
    selection = int(input("Select Project : ")) - 1
    os.chdir(projects[selection].name + "/Source")
    print('\n---------Executing Selected Project---------')
    os.system("python master.py")
    print('\n-------------Execution Complete-------------')

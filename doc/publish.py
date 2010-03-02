"""
Script to automate publishing of cmepy documentation to github
"""

import os
import os.path
import subprocess
import optparse

def git_branch_name():
    """
    Returns the current git branch name.
    """
    
    p = subprocess.Popen(
        ['git', 'branch', '--no-color', '--contains', 'HEAD'],
        stdout=subprocess.PIPE
    )
    p.wait()
    output, _ = p.communicate()

    try:
        branch_name = (output.strip().split('*')[1]).strip()
        return branch_name
    except IndexError:
        raise RuntimeError('unable to parse current branch name')

def git_uncommited_changes():
    """
    Returns True if the current git branch has uncommitted changes.
    """
    p = subprocess.Popen(
        ['git', 'status'],
    )
    returncode = p.wait()
    return (returncode is not 0)

def execute(args, cwd = None, must_return_zero = True):
    """
    run args as a command through the shell
    """
    command = ' '.join(args)
    p = subprocess.Popen(command, shell = True, cwd = cwd)
    return_code = p.wait()
    if must_return_zero and (return_code is not 0):
        raise RuntimeError('failure: return_code %d' % return_code)
    

def publish_documentation(publish_branch, build_dir, remote_repo):
    """
    Publishes documentation
    
    Builds the documentation, moves it to a new branch named publish_branch,
    then pushes it to the remote_repo
    """
    
    def announce(msg):
        print '*** [ '+msg+' ]'
    
    current_branch = git_branch_name()
    announce('current branch is: \'%s\'' % current_branch)
    
    current_dir = os.path.basename(os.getcwd())
    
    announce('current directory is: \'%s\'' % current_dir)
    
    if publish_branch == current_branch:
        raise ValueError('publish_branch cannot be current_branch, aborting')
    
    if git_uncommited_changes():
        raise RuntimeError('uncommitted changes in current_branch, aborting')
    
    announce('removing local publish_branch \'%s\'' % publish_branch)
    
    execute(['git', 'branch', '-D', publish_branch], must_return_zero = False)
    
    announce('building html documentation')
    
    execute(['make', 'clean'])
    execute(['make', 'html'])
    
    announce('exporting documentation to publish_branch')
    execute(
        [
            'git',
            'symbolic-ref',
            'HEAD',
            os.path.join('refs', 'heads', publish_branch)
        ],
        cwd = os.path.pardir
    )
    execute(
        [
            'rm',
            os.path.join('.git', 'index')
        ],
        cwd = os.path.pardir
    )
    
    build_path = os.path.join(current_dir, build_dir, 'html')
    # nb force the add with -f as build_path is in .gitignore file
    execute(
        ['git', 'add', '-f', build_path],
        cwd = os.path.pardir)
    execute(
        ['git', 'clean', '-fdx'],
        cwd = os.path.pardir
    )
    execute(
        ['git', 'mv', os.path.join(build_path, '*'), os.path.curdir],
        cwd = os.path.pardir
    )
    execute(
        ['rm', '-rf', current_dir],
        cwd = os.path.pardir
    )
    # add a .nojekyll file to stop github using jekyll, as jekyll fails to
    # render pages / directories with a leading underscore
    nojekyll = '.nojekyll'
    execute(
        ['touch', nojekyll],
        cwd = os.path.pardir
    )
    execute(
        ['git', 'add', nojekyll],
        cwd = os.path.pardir
    )
    announce('making a local commit')
    commit_msg = '\'automatically rebuilt documentation branch\''
    execute(
        ['git', 'commit', '-m', commit_msg],
        cwd = os.path.pardir
    )
    announce('removing publish_branch from remote_repo')
    execute(
        ['git', 'push', remote_repo, ':' + publish_branch],
        cwd = os.path.pardir,
        must_return_zero = False,
    )
    announce('pushing publish_branch to remote_repo')
    execute(
        ['git', 'push', remote_repo, publish_branch],
        cwd = os.path.pardir
    )
    announce('all finished!')
    
def main():
    """
    a nice command line interface for publish_documentation()
    """
    parser = optparse.OptionParser()
    parser.set_defaults(
        publish_branch = 'gh-pages',
        build_dir = '_build',
        remote_repo = 'origin'
    )
    parser.add_option(
        '-p',
        '--publish-branch',
        action='store',
        type='string',
        dest='publish_branch'
    )
    parser.add_option(
        '-b',
        '--build-dir',
        action='store',
        type='string',
        dest='build_dir'
    )
    parser.add_option(
        '-o',
        '--origin',
        action='store',
        type='string',
        dest='remote_repo'
    )
    
    (options, args) = parser.parse_args()
    if len(args) != 0:
        parser.error('incorrect number of arguments')
    
    publish_documentation(
        options.publish_branch,
        options.build_dir,
        options.remote_repo
    )
    
if __name__ == '__main__':
    main()

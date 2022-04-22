## The basics:

#### 1. Create  a local copy of the default branch of a repository.
```
git clone https://gitlab.gfdl.noaa.gov/fms/{repo_name}.git 
```
If you are running from an xml, you can skip this step and cd to your source directory. 

#### 2. Create a new local branch to make changes on
Your branch should be in the format of /user/**ini**/**BranchName** where **ini** 
are your initals and **BranchName** is a name that describes the purpose of your
branch.  Avoid using a BranchName like *bugfix* or *update* because this is not
descritive.
```
git checkout -b {branch name}
```
#### 3. Make code changes
#### 4. Commit changes

##### 4.1 See the changes that you made:
```
git status
```

##### 4.2 Add changes to the commit:
```
git add path/file
```
which will add only the files specified. 

```
git add -u 
```
which will add only the files that were changed and already exist in the repo.

```
git add .
```
which will add all the files that were changed including any new files.

##### 4.3 Create the commit:
```
git commit -m "Descriptive message of what was done"
```
#### 5. Push changes
```
git push origin branch_name
```

#### 6. Repeat steps 3-6 until code is finished!

## Creating a merge request:
Once your code updates at pushed to gitlab, you are ready to create a merge request!
#### 1. Pulling changes from main branch:
Before submitted a merge request, it is important to keep your development branch in synch with the main repo, this will help avoid any merge conflicts that may arise. 
```
git fetch
git pull origin main
```
#### 1.2 Fixing merge conflicts:
At some point, you will receive messages that files have merge conflicts after pulling changes from the main branch. This happens when both branches/commits changed the same lines of code. These must be resolved before your code is merged to the main branch. 
- `git status` will show files that need to be resolved. 
- Open the file with your favorite text editor
- Fix the conflicts and delete the <<<<< ====== >>>>>> lines
- Add your files `git add file`
- Commit your changes

#### 2. Opening a merge request:
[Submit a merge request using the web interface](https://gitlab.gfdl.noaa.gov/fms/am5_phys/-/merge_requests/new)

- The source branch is your branch.
- The target branch is the `main` branch.
- The merge request tile should be short and descriptive
- Add a description of what you changed
- If there is a issue that this merge request is solving link it to the merge request
- Indicate how the code was tested. What xml was used? What experiments did you used? What compiler? What system? And any other information.
- Indicate whether your answers reproduce and any namelist options that need to be added. If answers do not reproduce provide an explanation. 
- Indicate at least one reviewer to review your code
- Complete the checklist
- Once the reviewer(s) approves assign the merge request to Uriel Ramirez

## Reporting bugs or issues
[Open a new issue using the web interface](https://gitlab.gfdl.noaa.gov/fms/am5_phys/-/issues/new) describing the bug you are solving or the feature you are adding to the code. 

- The issue title should be short and descriptive. 
- The issue description should be clear and concise. Include enough information to help others reproduce the issue, or understand the change requested. 
- Assign the issue to the person that will be fixing it.

## Other helpful commands
```
git diff
```
This will show you the uncomitted changes.

```
git diff <commit hash> <commit hash> 
```
This will show you the difference between two commits

```
git revert <commit hash> 
```
This will undo a commit, while keeping past commits. It will also create a new commit showing that it has been reverted. Best for keeping the history intact when the reverted commits have already been pushed to the remote repository.

```
git reset <commit hash>
```
This will go back to the given commit while deleting any past commits after the given hash. This essentially sets the repository back to it's state from the given commit. Best for when commits only exist locally and can be lost without changing the remote's commit history.
```
git cherry-pick <commit hash>
```
This will apply an existing commit from another branch to your current branch, a commit can be 'cherrypicked' and added on top of the current history. Similar to a merge, but only adds one commit to the top instead of merging an entire history.

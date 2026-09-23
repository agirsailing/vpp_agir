# Velocity Prediction Program (VPP)

A velocity prediction program (VPP) to analyse the performance of the Moth.

The [hull hydrostatics guide](HYDROSTATICS.md) covers explicit `.bri` import,
upright and inclined hydrostatics, equilibrium draft, and free-floating attitude.
It includes an independently validated box example and MATLAB unit tests.

## Contributing

This guide explains how members of the Ägir sailing team can contribute to the VPP. Work on a separate branch for each feature or fix, then open a pull request so the team can review your changes.

### 1. Clone the repository

Clone the repository to create a local copy on your computer. You only need to do this once.

```bash
git clone https://github.com/agirsailing/vpp_agir.git
cd vpp_agir
```

### 2. Create a feature branch

Before starting new work, make sure your local `main` branch is up to date, then create a branch with a short, descriptive name:

```bash
git switch main
git pull --ff-only
git switch -c feature/add-drag-model
```

Use a name that describes your work, such as `feature/add-drag-model` or `fix/foil-lift-calculation`.

### 3. Make and commit your changes

Make your changes and check that the affected code works. Review your changes before committing:

```bash
git status
git diff
```

Stage the files you want to include:

```bash
git add path/to/file.py
```

or, to stage all file changes:

```bash
git add .
```

Then commit them with a clear message describing the change:

```bash
git commit -m "Add initial foil drag model"
```

Replace `path/to/file.py` with the actual file path. You can stage multiple files and make several commits on the same branch.

### 4. Push your branch

Push your branch to GitHub:

```bash
git push -u origin feature/add-drag-model
```

After the first push, you can upload further commits with `git push`.

### 5. Open a pull request

On [GitHub](https://github.com/agirsailing/vpp_agir), open a pull request from your feature branch into `main`.

Give it a clear title and briefly describe:
- What you changed and why.
- How you tested the changes.
- Anything reviewers should pay particular attention to.

Ask another team member to review your pull request before merging. If changes are requested, commit and push them to the same branch—the pull request will update automatically.

#!/bin/env node

const fs = require("fs");
const semver = require("semver");
const child_process = require("child_process");

function getTagVersionFromNpm(tag) {
    try {
        return child_process
            .execSync(`npm info ${package.name} version --tag="${tag}"`)
            .toString("utf8")
            .trim();
    } catch (e) {
        return null;
    }
}

function setVersion(version) {
    child_process.execSync(`npm version ${version} --no-git-tag-version --allow-same-version`);
}

// load package.json
const package = JSON.parse(fs.readFileSync("package.json", "utf8"));
console.log(`Package metadata read. Version: ${package.version}`);
// get the latest version from npm registry
const currentLatest = getTagVersionFromNpm("latest") || "0.0.0";
console.log(`Artifactor version: ${currentLatest}`);
// increment the next version number
const nextVersion = semver.inc(currentLatest, "prerelease", "ccdc");
console.log(`NextVersion version: ${nextVersion}`);

const nextPackageVersion = semver.inc(package.version, "prerelease", "ccdc");
// check if the package.json version is greater than the next version
const publishTag = semver.gt(nextPackageVersion, nextVersion, { includePrerelease: true })
    ? nextPackageVersion
    : nextVersion;
console.log(`PublishTag: ${publishTag}`);

// set the version in package.json
setVersion(publishTag);
